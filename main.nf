/*
 * Flujo de trabajo
 * En este archivo se ejecutan los metodos de busqueda de fusiones:
 * -RNApeg
 * -Cicero
 * -Arriba
 * -FusionCatcher
 * -Rascall
 * Cuantifica con Salmon.
 *
 */
nextflow.enable.dsl=2


include {
          STAR_aligner;
          RNApeg;
          Cicero;
          Salmon;
          Arriba;
          FusionCatcher;
          FusionList;
	        Fungi;
	        Rascall;
          FreeBayes;
          SnpEff;
          FusionSummary;
          FastQC as FastQC_before;
          FastQC as FastQC_after;
          Fastp;
          MultiQC as MultiQC_before;
          MultiQC as MultiQC_after;
          } from './modules/modules.nf'

workflow {

  
  /*
   * Genera el canal "fqs_ch" a partir del archivo ".tsv", 
   * contiene las columnas "Sample", "R1" y "R2".  
  */
  fqs_ch = Channel.fromPath(file(params.runSampleSheet))
  				.splitCsv(header: true, sep: '\t')
  				.map { sample -> [sample["Sample"], file(sample["R1"]), file(sample["R2"])]}
  
  // Control de calidad antes del filtrado
  FastQC_before(fqs_ch, "beforeTrimm")

  // MultiQC para lecturas sin filtrar
  MultiQC_before(FastQC_before.out.qc.collect(), params.reportsDir+"/beforeTrimm")

  // Eliminacion de adaptadores y filtrado de calidad
  Fastp(fqs_ch)

  // Control de calidad despues del filtrado
  FastQC_after(Fastp.out.reads, "afterTrimm")

  // Salmon cuantifica las muestras.
  Salmon(params.referenceDir, Fastp.out.reads)

  // a partir de Fastp.out.reads, Rascall busca fusiones.
  Rascall(Fastp.out.reads)

  // STAR se utiliza en su modo de alineamiento.
  STAR_aligner(params.referenceDir, params.threadsSTAR, Fastp.out.reads)

  // FreeBayes llama variantes a partir del BAM de STAR.
  FreeBayes(params.referenceDir, STAR_aligner.out.BAM)

  // SnpEff anota las variantes con informacion funcional.
  SnpEff(FreeBayes.out.vcf)

  // RNApeg se utiliza para busqueda de uniones y fusiones, utiliza los ".bam" resultado de STAR_aligner.
  RNApeg(params.referenceDir,STAR_aligner.out.BAM)
  
  // Se genera un canal para cicero, este tiene la union de los ".bam" de STAR y las uniones de RNApeg.
  cicero_ch = STAR_aligner.out.BAM
    .join(RNApeg.out.junctions)
  
  // Cicero busca fusiones a partir del canal previamente creado y genera un reporte
  Cicero(params.referenceDir,cicero_ch)
  
  // Arriba busca fusiones utilizando los ".bam" resultado de STAR_aligner, genera un reporte con las fusiones.
  Arriba(params.referenceDir,STAR_aligner.out.BAM)

  // FusionCatcher busca fusiones y genera un reporte.
  FusionCatcher(params.referenceDir, Fastp.out.reads)

  /* Genera un canal con los resultados de los reportes de fusiones de FusionCatcher, 
   * Arriba y Cicero. No avanza hasta que se obtienen todos los reportes.
   */
  fusions_ch = Cicero.out.fusions
                     .concat(Arriba.out.fusions, FusionCatcher.out.fusions)
                     .collect()

  // FusionList genera una tabla con los resultados por muestra de los metodos de busqueda de Fusiones (FusionCatcher, Arriba y Cicero) para que Fungi los pueda procesar.
  FusionList(fusions_ch)  

  // Fungi analiza y genera un consenso de la lista de FusionList.
  Fungi(FusionList.out.list)
   
  // FusionSummary genera un reporte usando los reportes de Cicero, Arriba, FusionCatcher y Rascall.
  FusionSummary(params.runSampleSheet, Fungi.out.bp, Rascall.out.results.collect())

  // MultiQC final: FastQC (antes y despues) + Fastp + STAR
  reports_ch = FastQC_before.out.qc
                 .mix(FastQC_after.out.qc)
                 .mix(Fastp.out.json)
                 .mix(STAR_aligner.out.logs)
                 .collect()

  MultiQC_after(reports_ch, params.reportsDir+"/afterTrimm")

}
