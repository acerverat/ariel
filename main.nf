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
 * Expresion Cluster genera clusters, matriz de TPMs y log2TPM
 */
nextflow.enable.dsl=2


include { 
          STAR_aligner; 
          RNApeg;
          Cicero;
          Salmon;
          ExprClusters;
          Arriba;
          FusionCatcher;
          FusionList;
	        Fungi;
	        Rascall;
          FusionSummary;
	      //  PreExprCluster;
          } from './modules/modules.nf'

workflow {

  
  /*
   * Genera el canal "fqs_ch" a partir del archivo ".tsv", 
   * contiene las columnas "Sample", "R1" y "R2".  
  */
  fqs_ch = Channel.fromPath(file(params.runSampleSheet))
  				.splitCsv(header: true, sep: '\t')
  				.map { sample -> [sample["Sample"], file(sample["R1"]), file(sample["R2"])]}
  
   // First quality check out
  FastQC(fqs_ch, "beforeTrimm")


  // MultiQC for unfiltered reads
  MultiQC(FastQC.out.qc.collect(), params.reportsDir+"/beforeTrimm")

   // Adaptors elimination and qualiy filtering
  Fastp(fqs_ch)

  // Quality check out after trimming
  FastQC(Fastp.out.reads, "afterTrimm")

  // Salmon cuantifica las muestras.
  Salmon(params.referenceDir,fqs_ch)

  // "salmon_ch" espera a obtener todos los resultados antes que "PreExprCluster" los reciba.
  salmon_ch = Salmon.out.sf.collect()

  // ExprClusters fenera la matriz de expresion
  ExprClusters(params.runSampleSheet, salmon_ch, params.ensg_enst_table, params.exprDir)
  
  // a partir de fqs_ch, Rascall busca fusiones.
  Rascall(fqs_ch)

  // STAR se utiliza en su modo de alineamiento. 
  STAR_aligner(params.referenceDir,params.threadsSTAR,fqs_ch) 

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
  FusionCatcher(params.referenceDir,fqs_ch)

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
   
  Rascall.out.results.view()
  ExprClusters.out.clusters.view()  
  // FusionSummary genera un reporte usando los reportes de Cicero, Arriba, FusionCatcher y Rascall. 
  FusionSummary(params.runSampleSheet,Fungi.out.bp,Rascall.out.results.collect(),ExprClusters.out.clusters,params.method_counts,params.supporting_reads)

    // MultiQC for Fastp+Fastp+Salmon+STAR
  reports_ch = FastQC.out.qc.collect().concat(Fastp.out.reads.collect())
                                      .concat(Salmon.out.quant.collect())
                                      .concat(STAR_aligner.out.BAM.collect()).collect()

  MultiQC(reports_ch, params.reportsDir+"/afterTrimm")

}
