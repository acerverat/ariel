  /*
   *                        		---- MODULOS ----
   *
   * En este archivo se encuentran los modulos que se utilizan en el flujo de trabajo.
   * En cada uno se especifica su input, ouput y una descripcion.
   * 
   * 
   *       
   */


process FusionSummary {
  /*
   *                        ---- FusionSummary ----
   *
   * FusionSummary genera los reportes finales del pipeline en dos pasos secuenciales.
   *
   * Input:
   *   - SampleSheet (path): TSV con las muestras a analizar.
   *   - bp_consensus (path): Reporte de Fungi con los breakpoints.
   *   - rascall (val): Archivos de resultados de RaScALL (dependencia de ejecucion).
   *   - cluster (path): Tabla log2tpm_CRLF2_3_clusters.tsv de ExprClusters.
   *
   * Output:
   *   - hallazgos_principales.csv:
   *       Muestra, Fusion, Metodos, Subtipo, Subtipo_Emergente, Punto_de_corte,
   *       SR_Arriba, SR_Cicero, SR_Fusioncatcher
   *
   *   - hallazgos_otros.csv:
   *       Muestra, Fusion, Metodos, Subtipo, Subtipo_Emergente, CRLF2_expr,
   *       SNV_RaScALL, Deleciones_Focales_Rascall, DUX4r_Rascall, Duplicacion_CRLF2
   */
  cache 'lenient'
  publishDir params.resultsDir+"/reports", mode: 'copy'

  input:
    path SampleSheet
    path bp_consensus
    val rascall
    path cluster

  output:
    file ('hallazgos_*.csv')

  script:
  """
    export TMPDIR=\$PWD

    # Paso 1: genera hallazgos_principales.csv, fusiones_otras.csv y rascall_data.csv
    generaReporteHallazgosPrincipales.R ${SampleSheet} ${bp_consensus} ${params.resultsDir}/rascall

    # Paso 2: genera hallazgos_otros.csv usando los intermedios del paso anterior
    generaReporteHallazgosOtros.R fusiones_otras.csv rascall_data.csv ${params.resultsDir}/fusions/cicero ${SampleSheet} ${cluster}
  """
}

process ExprClusters {
  /*
   *                        ---- ExprClusters ----
   * 
   * ExprClusters agrupa muestras en clusters de acuerdo con su expresion y genera graficos
   * para su visualizacion. Utiliza el contenedor de Docker 'ariel-env'.
   *
   *  Input:
   *	- runSampleSheet (path): Directorio con las rutas de las muestras.
   *	- SalmonCollect (val): Lista con las rutas de los archivos de cuantificacion de Salmon.
   *	- tablaGenesReferencia (path): Tabla con el gene id y ensembl id para poder unir la tabla generada de Salmon collect y tablaTPM.
   *    - tablaTPM (path): Tabla de TPMs de muestras de un panel de genes.   
   *
   *  Output:
   *  	- log2tpm_CRLF2_3_clusters.tsv (path): matriz transformada a log(TPM + 1)
   *	- log2tpm_{gen}_{k}_clusters.tsv (path): tabla con: Sample, expresion y cluster
   *	- boxplot_{gen}.png (path): boxplot con expresion por cluster.
   *	- boxplot_nombres_{gen}: boxplot con nombre de las muestras seleccionadas. 
   *  
   */
  cache 'lenient'
  publishDir params.resultsDir+"/ExprClusters", mode: 'copy'
  publishDir params.resultsDir+"/reports", mode: 'copy', pattern: "*png"
  publishDir params.resultsDir+"/reports", mode: 'copy', pattern: "*cluster*"
  container 'ariel-env:latest'

  input:
    path runSampleSheet
    val SalmonCollect
    path tablaGenesReferencia
    path tablaTPM	

  output:
    path("log2tpm_CRLF2_3_clusters.tsv"), emit: clusters
    path("*tpm.tsv")
    path("*png")

  script:
  """
    export TMPDIR=\$PWD

    # filtra caracteres para que preExprCluster.R pueda leer las rutas
    salmoncsv=\$(echo ${SalmonCollect} | tr -d ' ' | tr -d '[' | tr -d ']')

    # ejecuta preExprCluster para obtener la tabla de tpms
    preExprCluster.R \${salmoncsv} ${tablaGenesReferencia} ${tablaTPM}    

    # ejecuta kmeans con los genes CRLF2 y DUX4
    kmeans.R ${runSampleSheet} tpm.tsv  CRLF2 3
    kmeans.R ${runSampleSheet} tpm.tsv  DUX4 3
  """
}

process Rascall {
 /*
  *                        ---- Rascall ----
  * 
  * RaScALL detecta fusiones, SNVs, fusiones IGH, deleciones focales y DUX4.
  * Comparado con los otros metodos de busqueda de fusiones, este se especializa
  * en Leucemia Linfoblastica Aguda, usando su propia base de datos.
  * Rascall utiliza el contenedor de Docker: rascall:1.0
  * Genera un reporte general y varios reportes en su directorio de trabajo.
  *
  * Input:
  *   - tuple:
  *       - sample (val): Nombre de la muestra.
  *       - R1 (file): Lecturas forward (FASTQ).
  *       - R2 (file): Lecturas reverse (FASTQ).
  *
  * Output:
  *   - tuple (emit: results):
  *       - sample (val): Nombre de la muestra.
  *       - rascall_dir/*//*final_variants.csv (file): 
  *			-Reporte con los campos:
  *				File, Target_Type, Alteration, Query, Type, Variant_name,
  *				rVAF, Expression, Min_coverage, Sequence, Reference_sequence.
  *
  */
  cache 'lenient'
  publishDir params.resultsDir+"/rascall", mode: 'copy'

  input:
  tuple val(sample), file(R1), file(R2)

  output:
  tuple val(sample), file("rascall_dir/*/*final_variants.csv"), emit: results

  script:
  """
    # se genera directorio para RaScALL
    mkdir -p rascall_dir/${sample}

    # directorio real donde se encuentran los FASTQ
    data=\$( dirname \$( readlink -f ${R1} ) )

    # nombres de los archivos sin ruta
    r1=\$( basename ${R1} )
    r2=\$( basename ${R2} )

    # se utiliza el contenedor de docker que contiene RaScALL
    # WORKDIR del contenedor es /opt/rascall, por lo que run_km.sh resuelve:
    #   VIRTUAL_ENV=/opt/rascall/.virtualenvs/km  y  OUTPUTS=/opt/rascall/output
    # rascall_dir/ del host se monta en /opt/rascall/output/ para que los resultados
    # queden disponibles en el directorio de trabajo de Nextflow
    docker run -u \$(id -u):\$(id -g) --rm \
               -w /opt/rascall \
               -e HOME=/tmp \
               -v \${data}:/data \
               -v \$PWD/rascall_dir:/opt/rascall/output \
               rascall:1.0 \
    -c "ln -sf /data/\${r1} /opt/rascall/output/${sample}/${sample}_R1.fastq.gz;
        ln -sf /data/\${r2} /opt/rascall/output/${sample}/${sample}_R2.fastq.gz;
        bash /RaScALL/run_km.sh /opt/rascall/output/${sample}/${sample}_R1.fastq.gz /opt/rascall/output/${sample}/${sample}_R2.fastq.gz ${params.threadsRascall} /RaScALL/ALL_targets/DUX4;
        bash /RaScALL/run_km.sh /opt/rascall/output/${sample}/${sample}_R1.fastq.gz /opt/rascall/output/${sample}/${sample}_R2.fastq.gz ${params.threadsRascall} /RaScALL/ALL_targets/Fusion;
        bash /RaScALL/run_km.sh /opt/rascall/output/${sample}/${sample}_R1.fastq.gz /opt/rascall/output/${sample}/${sample}_R2.fastq.gz ${params.threadsRascall} /RaScALL/ALL_targets/IGH_fusion;
        bash /RaScALL/run_km.sh /opt/rascall/output/${sample}/${sample}_R1.fastq.gz /opt/rascall/output/${sample}/${sample}_R2.fastq.gz ${params.threadsRascall} /RaScALL/ALL_targets/SNV;
        bash /RaScALL/run_km.sh /opt/rascall/output/${sample}/${sample}_R1.fastq.gz /opt/rascall/output/${sample}/${sample}_R2.fastq.gz ${params.threadsRascall} /RaScALL/ALL_targets/focal_deletions;
        Rscript /RaScALL/bin/filter_km_output.R /opt/rascall/output/${sample}"
  """
}

process Fungi { 
  
 /*
  *                        ---- Fungi ----
  * 
  * Fungi utiliza scripts del repositorio de Fungi de la Dra. Alejandra Cervera.
  * Primero analiza los reportes de Arriba, FusionCatcher y Cicero, y despues hace un
  * consenso.
  * Utiliza los scripts fungi-fusion-analyzer y fungi-fusion-consensus de Fungi instalado
  * en el contenedor de Docker ariel-env.
  * 
  *
  * Input:
  *   - tuple:
  *       - sample (val): Nombre de la muestra.
  *       - R1 (file): Lecturas forward (FASTQ).
  *       - R2 (file): Lecturas reverse (FASTQ).
  *
  * Output:
  *   - tuple (emit: results):
  *       - bp_consensus_report.tsv (path): 
  *		Reporte de los breakpoints con los campos:
  *			FusionName, Sample, best_bp Methods_count, Methods_list, Supporting_reads, same_score_bp, 
  * 			Sample_occurrence_method_score, SampleCount, best_by_sample_count, same_score_sample_count,
  *			Sample_occurence_sample_score, annotations, additional_info
  *
  *       - combined_fusions_report.tsv (path):
  *		Reporte con los campos:
  *			FusionName, bp_coordinates, Sample, Supporting_reads, Methods Sample_occurrence, annotations,
  *			additional_info, MethodsCount, SampleCount
  *
  *
  */
  cache 'lenient' 
  publishDir params.resultsDir+"/fungi", mode: 'copy'

  input:
    val list

  output:
    path("fungi_output/consensus/bp_consensus_report.tsv"), emit: bp
    path("fungi_output/consensus/combined_fusions_report.tsv")
  script:
    """
    outdir=\$PWD
    fcdb="${params.referenceDir}/fusioncatcher_db"

    # analisis de fusiones
    # se monta fusioncatcher_db sobre /opt/fusioncatcher/data para que exons.txt
    # (generado durante la descarga de la db) este disponible para fungi
    docker run -t --rm \
    -v ${workDir}:${workDir} \
    -v \${fcdb}:/opt/fusioncatcher/data \
    ariel-env:latest \
    fungi-fusion-analyzer -c myConfig.txt -o \${outdir}/fungi_output --input-list ${list} annotate --filter-ensembl 'invalid_gene,same_gene,homologs' --filter-db 'banned,paralog' --filter-min-count 0

    # consenso
    docker run -t --rm \
    -v ${workDir}:${workDir} \
    -v \${fcdb}:/opt/fusioncatcher/data \
    ariel-env:latest \
    fungi-fusion-consensus -o \${outdir}/fungi_output/consensus --fungi_annotated \${outdir}/fungi_output/annotated
    """
}


process FusionList {
/*
 *                        ---- FusionList ----
 * 
 * Convierte una lista de fusiones en formato CSV de los metodos de busqueda de fusiones
 * Arriba, FusionCatcher y Cicero a un TSV compatible con Fungi.
 *
 * Input:
 *   - fusionList (val): Lista plana con elementos en grupos de 3:
 *       [sample, tool, file, sample, tool, file, ...]
 *
 * Output:
 *   - input-list.txt (path):
 *   Archivo TSV con formato:
 *         CA001	cicero	/ruta...
 *         CA002  arriba   /ruta...
 *         CA003  fusioncatcher  /ruta...
 *
 */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir+"/fusions", mode: 'copy'

  input:
    val fusionList

  output:
    path 'input-list.txt', emit: list

  script:                    
  """
    # Convierte fusionList de una lista con las fusiones por metodo a un formato tsv para que Fungi lo pueda recibir
    echo -e "Sample\tTool\tFile" > input-list.txt
    echo "${fusionList.join(',')}" | sed -E 's/(,[^,]*,[^,]*),/\\1\\n/g' | tr ',' '\\t' >> input-list.txt
  """
}

process Salmon {
  /*
   *                        ---- Salmon ----
   * 
   * Salmon realiza la cuantificacion de expresion de genes. Utiliza la version
   * de Salmon instalada en el contenedor de Docker 'ariel-env'.
   *
   * Input:
   *   - referenceDir (path): Directorio donde esta el indice generado por Salmon_index.
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - R1 (file): Lecturas forward (FASTQ).
   *       - R2 (file): Lecturas reverse (FASTQ).
   *
   * Output:
   *   - {sample}.quant.sf (path): 
   *      Archivo de cuantificacion en formato '.sf' separado por tabs con las siguientes columnas: 
   *      - Name  Length  EffectiveLength TPM     NumReads
   *
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir+"/quantification", mode: 'copy'

  input:
    path referenceDir
    tuple val(sample), file(R1), file(R2)

  output:
    path '*.sf', emit: sf

  script:
  """
    # se utiliza el indice de Salmon ubicado en el directorio de referencia
    transcriptome_idx=\$PWD/${referenceDir}/salmon_index
    
    # Ejecuta Salmon para cuantificar 
    salmon quant -i \${transcriptome_idx} -l A \
                 -1 ${R1} \
                 -2 ${R2} \
                 -p 4 --validateMappings -o \$PWD/${sample}_quant
              
      # copia los resultados al directorio de trabajo
      cp \$PWD/${sample}_quant/quant.sf \$PWD/${sample}.quant.sf
  """
}

process STAR_aligner {
  /*
   *                        ---- STAR_aligner ----
   *
   * STAR_aligner realiza el alineamiento de las lecturs de RNA contra un genoma de
   * referencia, generando archivos BAM y BAI requeridos por Cicero, RNApeg y Arriba.
   *
   * Input:
   *   - referenceDir (path): Directorio con el indice generado por STAR_index.
   *   - threadsSTAR (val): Numero de hilos que usara STAR_aligner.
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - R1 (file): Lecturas forward (FASTQ).
   *       - R2 (file): Lecturas reverse (FASTQ).
   *
   * Output:
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - {sample}_Aligned.sortedByCoord.out.bam (file): Archivo BAM de la muestra con el genoma de referencia.
   *       - {sample}_Aligned.sortedByCoord.out.bam.bai (file): Archivo BAI indice del BAM. 
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir+"/alignments", mode: 'copy'
  beforeScript 'chmod o+rw .'

  input:
    path referenceDir
    val threadsSTAR
    tuple val(sample), file(R1), file(R2)
  
  output:
    tuple val(sample),
          file("${sample}_Aligned.sortedByCoord.out.bam"),
          file("${sample}_Aligned.sortedByCoord.out.bam.bai"),
          emit: BAM
    path("${sample}_Log.final.out"), emit: logs


  script:
  """
    set -eou pipefail

    # directorio de referencia con el indice de STAR
    genome_idx=\$PWD/${referenceDir}/STAR_2.7.10b_index

    # Corre STAR en modo alineamiento
    STAR --runMode alignReads \
       --genomeDir  \${genome_idx} \
       --runThreadN ${threadsSTAR} \
       --readFilesIn ${R1} ${R2} \
       --outFileNamePrefix ${sample}"_" \
       --outReadsUnmapped None \
       --twopassMode Basic \
       --twopass1readsN -1 \
       --readFilesCommand "gunzip -c" \
       --outSAMunmapped Within \
       --outSAMtype BAM SortedByCoordinate \
       --outBAMcompression 0 \
       --outTmpDir ./tmp_star_${sample} \
       --outSAMattributes NH HI NM MD AS nM jM jI XS \
       --chimSegmentMin 20 \
       --chimJunctionOverhangMin 20 \
       --chimOutType WithinBAM

     samtools index ${sample}_Aligned.sortedByCoord.out.bam
  """
}



process Arriba {
  /*
   *                        ---- Arriba ----
   * ARRIBA busca fusiones y genera un reporte. Ademas, produce un pdf con
   * los dibujos de las fusiones encontradas.
   * Utiliza la version de Arriba descargada en el contenedor 'ariel-env'
   *
   * Input:
   *  - referenceDir (path): Directorio con el genoma de referencia que utiliza Cicero y sus anotaciones
   *    - se pueden obtener desde el repositorio de cicero:
   *         https://github.com/stjude/Cicero?tab=readme-ov-file#downloading-reference-files-
   *
   *  - tuple:
   *    - sample (val): Nombre de la muestra.
   *    - bamfile (file): Archivo BAM.
   *    - bai (file): Indice del BAM.
   *
   * Output:
   *   - tupla:
   *       - sample (val) : Nombre de la muestra-
   *       - 'arriba' (val) : Etiqueta del metodo para el modulo Fungi-
   *       - file: "${sample}_fusions.tsv" Reporte de fusiones encontradas.
   *
   *   - {sample}_fusions.pdf:  Dibujos de las fusiones encontradas.
   *
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir + "/fusions/arriba", mode: 'copy'

  input:
    path referenceDir
    tuple val(sample), file(bamfile), file(baifile)

  output:
    tuple val(sample), val('arriba'), file("${sample}_fusions.tsv"), emit: fusions
    path("${sample}_fusions.pdf")

  script:
  """
    export TMPDIR=\$PWD

    # Agrega variables con los nombres de las rutas que requiere
    gtf=\$PWD/${referenceDir}/gencode.v42.annotation.gtf
    fasta=\$PWD/${referenceDir}/cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa
    arriba_home="/usr/local/lib/arriba_v2.4.0"

    # Correr Arriba
      "\${arriba_home}"/arriba \
      -x ${bamfile} \
      -g "\${gtf}" \
      -a "\${fasta}" \
      -b "\${arriba_home}"/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz \
      -k "\${arriba_home}"/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
      -p "\${arriba_home}"/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
      -o fusions.tsv

    # Dibujar Fusiones con Arriba
    "\${arriba_home}"/draw_fusions.R \
    --fusions=\$PWD/fusions.tsv \
    --alignments=${bamfile} \
    --output=fusions.pdf \
    --annotation="\${gtf}"

    # Renombrar las fusiones con en nombre de la muestra
    mv \$PWD/fusions.pdf \$PWD/${sample}_fusions.pdf
    mv \$PWD/fusions.tsv \$PWD/${sample}_fusions.tsv

  """
}


process RNApeg {
  /*
   *                        ---- RNApeg ----
   * RNApeg genera el archivo de splice junctions que requiere Cicero.
   * Utiliza la ultima version de RNApeg.
   *
   * Input:
   *  - referenceDir (path): Directorio con el genoma de referencia y sus anotaciones
   *    - se pueden obtener desde el repositorio de RNApeg:
   *         https://github.com/stjude/RNApeg?tab=readme-ov-file#downloading-reference-files
   *    - tupla:
   *      - sample (val): Nombre de la muestra.
   *      - bamfile (file): Archivo BAM.
   *      - bai (file): Indice del BAM.
   *
   * Output:
   *   - tuple:
   *       - sample (val) : Nombre de la muestra.
   *       - file: "{bamfile}.junctions.tab.shifted.tab" Archivo de junctions.
   *
   */

  cache 'lenient'
  publishDir params.resultsDir+"/junction", mode: 'copy'
  input:
    path(referenceDir)
    tuple val(sample), file(bamfile), file(bai)

  output:
    tuple val(sample), file("${bamfile}.junctions.tab.shifted.tab"), emit: junctions

  script:
  """
    # obtiene el nombre, ruta y directorio del bamfile
    bam=\$( basename ${bamfile} )
    bamfile=\$( readlink ${bamfile} )
    bamdir=\$( dirname \${bamfile} )

    # crea las variables con las rutas del directorio de referencia y el directorio de trabajo
    refdir=\$( readlink ${referenceDir} )
    outdir=\$PWD

    # Ejecuta el contenedor de Docker, utilizando la ultima version de rnapeg
    docker run -u \$(id -u):\$(id -g) --rm \
                --mount type=bind,source=\${outdir},target=/results \
                -v "\${refdir}":/references \
                -v "\${bamdir}":/data \
                ghcr.io/stjude/rnapeg \
                    -b /data/\${bam} \
                    -f /references/cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa \
                    -r /references/cicero_references/Homo_sapiens/GRCh38_no_alt/mRNA/RefSeq/refFlat.txt
  """
}

process Cicero {
  /*
   *                        ---- Cicero ----
   *
   * CICERO detecta de fusiones a partir de un BAM, BAI y junctions, 
   * utilizando el contenedor oficial de Cicero.
   *
   * Input:
   *   - referenceDir (path): Directorio con el genoma humano de referencia.
   *      - las referencias se pueden obtener desde el repositorio de cicero:
   *        https://github.com/stjude/Cicero?tab=readme-ov-file#downloading-reference-files-
   *   - tupla:
   *       - sample (val): Nombre de la muestra.
   *       - bamfile (file): Archivo BAM.
   *       - bai (file): Indice del BAM.
   *       - junctions (file): Archivo junctions.
   *
   *  Output:
   *   - tupla:
   *       - sample (val): Nombre de la muestra.
   *       - 'cicero' (val): Etiqueta del metodo para Fungi.
   *       - "{sample}_final_fusions.txt" (file): Reporte de fusiones encontradas.
   *  
   *   - "${sample}_HQ_only_fusions.txt": Reporte fusiones de tipo medal o HQ.
   *
   *   - "${sample}_all_ITD.txt": Reporte fusiones tipo ITD o medal.
   *
   */
  cache 'lenient'
  publishDir params.resultsDir+"/fusions/cicero", mode: 'copy'

  input:
    path(referenceDir)
    tuple val(sample), file(bamfile), file(bai), file(junctions)

  output:
    tuple val(sample), val('cicero'), file("${sample}_final_fusions.txt"), emit:fusions
    path("${sample}_HQ_only_fusions.txt")
    path("${sample}_all_ITD.txt")

  script:
  """
    # input se genera para almacenar los bam, bai y junctions
    mkdir inputs

    # results se utiliza para guardar los resultados de cicero
    mkdir results

    # los archivos bam, bai y junctions se copian a /inputs
    bam=\$( basename ${bamfile} )
    cp -L \${bam} inputs/
    cp -L ${bai} inputs/

    junctions=\$( basename ${junctions} )
    cp -L \${junctions} inputs/

    # se generan variables con las rutas del directorio de referencia y el directorio de trabajo
    refdir=\$( readlink ${referenceDir} )
    outdir=\$PWD

    # Ejecuta docker, utilizando cicero desde su ultima version
    docker run -u \$(id -u):\$(id -g) --rm \
                --mount type=bind,source=\${outdir},target=/cicero \
                -v "\${refdir}":/references \
                ghcr.io/stjude/cicero:latest Cicero.sh \
                       -b /cicero/inputs/\${bam} \
                       -r /references/cicero_references \
                       -g GRCh38_no_alt \
                       -o /cicero/results/${sample} \
                       -j /cicero/inputs/\${junctions} \
                       -no-optimize -n 2

    # copia los resultados al directorio de trabajo
    cp results/${sample}/CICERO_DATADIR/*/final_fusions.txt \$PWD/${sample}_final_fusions.txt

    # filtra fusiones de calidad HQ o medal, excluyendo ITDs
    grep -E 'HQ|medal' \$PWD/${sample}_final_fusions.txt | grep -v ITD > \$PWD/${sample}_HQ_only_fusions.txt
    
    # Filtrar fusiones medal o de tipo ITDs
    grep -E 'medal|ITD' \$PWD/${sample}_final_fusions.txt > \$PWD/${sample}_all_ITD.txt
  """
}


process FusionCatcher {
  /*
   *                        ---- FusionCatcher ----
   * FusionCatcher se utiliza para la deteccion de fusiones utilizando
   * su base de datos.
   *
   * Input:
   *   - referenceDir (path): Directorio de referencias; debe contener fusioncatcher_db/
   *       con la base de datos descargada por generaReferencias.sh.
   *   - tuple:
   *       - sample (val): nombre de la muestra.
   *       - R1 (file): Lecturas forward (FASTQ).
   *       - R2 (file): Lecturas reverse (FASTQ).
   *
   * Output:
   *   - tuple (emit: fusions):
   *       - sample (val)
   *       - 'fusioncatcher' (val): Etiqueta del metodo para Fungi.
   *       - file: "${sample}_final-list_candidate-fusion-genes.txt"
   *
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir+"/fusions/fusioncatcher", mode: 'copy'

  input:
    path referenceDir
    tuple val(sample), file(R1), file(R2)

  output:
    tuple val(sample), val('fusioncatcher'), file("${sample}_final-list_candidate-fusion-genes.txt"), emit:fusions

  script:
  """
    # Crea rutas /input y /ouput en el directorio de trabajo
    mkdir \$PWD/input
    mkdir \$PWD/${sample}_output

    # Crea una ruta simbolica de R1 y R2 a input
    ln -s \$PWD/${R1} \$PWD/input/r1.fq.gz
    ln -s \$PWD/${R2} \$PWD/input/r2.fq.gz

    # Ejecuta FusionCatcher apuntando a la base de datos del directorio de referencias
    /opt/fusioncatcher/bin/fusioncatcher.py \
      -d \$PWD/${referenceDir}/fusioncatcher_db/current \
      -i \$PWD/input \
      -o \$PWD/${sample}_output

    # Los resultados se guardan con el nombre de la muestra al inicio
    cp \$PWD/${sample}_output/final-list_candidate-fusion-genes.txt \$PWD/${sample}_final-list_candidate-fusion-genes.txt
  """
}

process FastQC {
  /*
   *                        ---- FastQC ----
   *
   * FastQC realiza el control de calidad de las lecturas FASTQ, generando
   * reportes HTML y archivos ZIP procesables por MultiQC.
   *
   * Input:
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - R1 (file): Lecturas forward (FASTQ).
   *       - R2 (file): Lecturas reverse (FASTQ).
   *   - label (val): Punto de control ("beforeTrimm" o "afterTrimm").
   *
   * Output:
   *   - *_fastqc.zip (emit: qc): Archivos ZIP con reportes para MultiQC.
   *   - *_fastqc.html: Reportes HTML de FastQC.
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.reportsDir, mode: 'copy'

  input:
    tuple val(sample), file(R1), file(R2)
    val label

  output:
    path("*_fastqc.zip"), emit: qc
    path("*_fastqc.html")

  script:
  """
    fastqc --dir . ${R1} ${R2}
  """
}

process Fastp {
  /*
   *                        ---- Fastp ----
   *
   * Fastp elimina adaptadores y filtra lecturas por calidad.
   *
   * Input:
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - R1 (file): Lecturas forward (FASTQ).
   *       - R2 (file): Lecturas reverse (FASTQ).
   *
   * Output:
   *   - tuple (emit: reads):
   *       - sample (val): Nombre de la muestra.
   *       - R1_{sample}_fastp.fq.gz (file): Lecturas forward filtradas.
   *       - R2_{sample}_fastp.fq.gz (file): Lecturas reverse filtradas.
   *   - {sample}_fastp.json (emit: json): Reporte JSON para MultiQC.
   *   - {sample}_fastp.html: Reporte HTML de Fastp.
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir+"/trimmed", mode: 'copy', pattern: "*.fq.gz"
  publishDir params.reportsDir+"/afterTrimm", mode: 'copy', pattern: "*.{html,json}"

  input:
    tuple val(sample), file(R1), file(R2)

  output:
    tuple val(sample), file("R1_${sample}_fastp.fq.gz"), file("R2_${sample}_fastp.fq.gz"), emit: reads
    path("${sample}_fastp.json"), emit: json
    path("${sample}_fastp.html")

  script:
  """
    fastp \
      -i ${R1} -I ${R2} \
      -o R1_${sample}_fastp.fq.gz -O R2_${sample}_fastp.fq.gz \
      -h ${sample}_fastp.html -j ${sample}_fastp.json
  """
}

process FreeBayes {
  /*
   *                        ---- FreeBayes ----
   *
   * FreeBayes realiza el llamado de variantes (SNVs e indels) a partir del
   * BAM generado por STAR_aligner.
   *
   * Input:
   *   - referenceDir (path): Directorio de referencias; debe contener la FASTA
   *       en cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - bam (file): Archivo BAM alineado y ordenado.
   *       - bai (file): Indice del BAM.
   *
   * Output:
   *   - tuple (emit: vcf):
   *       - sample (val): Nombre de la muestra.
   *       - {sample}.vcf (path): Variantes en formato VCF.
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir params.resultsDir + "/variants/freebayes", mode: 'copy'

  input:
    path referenceDir
    tuple val(sample), file(bam), file(bai)

  output:
    tuple val(sample), path("${sample}.vcf"), emit: vcf

  script:
  """
    fasta=\$PWD/${referenceDir}/cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa

    freebayes \
      -f \${fasta} \
      --min-alternate-count 3 \
      --min-alternate-fraction 0.05 \
      --skip-coverage-above 5000 \
      ${bam} > ${sample}.vcf
  """
}

process SnpEff {
  /*
   *                        ---- SnpEff ----
   *
   * SnpEff anota las variantes generadas por FreeBayes con informacion
   * funcional (gen, efecto, impacto).
   * La base de datos debe descargarse previamente con scripts/references/07_snpeff_db.sh
   * y se monta en /snpeff_data dentro del contenedor.
   *
   * Input:
   *   - tuple:
   *       - sample (val): Nombre de la muestra.
   *       - vcf (path): Archivo VCF de FreeBayes.
   *
   * Output:
   *   - tuple (emit: vcf):
   *       - sample (val): Nombre de la muestra.
   *       - {sample}_annotated.vcf (path): VCF anotado.
   *   - {sample}_snpeff_summary.html: Reporte HTML de SnpEff.
   *   - {sample}_snpeff_genes.txt: Tabla de genes anotados.
   */
  cache 'lenient'
  container 'ariel-env:latest'
  containerOptions "-v ${file(params.referenceDir).toAbsolutePath()}/snpeff_db:/snpeff_data"
  publishDir params.resultsDir + "/variants/snpeff", mode: 'copy'

  input:
    tuple val(sample), path(vcf)

  output:
    tuple val(sample), path("${sample}_annotated.vcf"), emit: vcf
    path("${sample}_snpeff_summary.html")
    path("${sample}_snpeff_genes.txt")

  script:
  """
    export TMPDIR=\$PWD

    java -Xmx4g -jar /opt/snpEff/snpEff.jar \
      -dataDir /snpeff_data \
      -stats ${sample}_snpeff_summary.html \
      ${params.snpeffGenome} \
      ${vcf} > ${sample}_annotated.vcf

    mv snpEff_genes.txt ${sample}_snpeff_genes.txt
  """
}

process MultiQC {
  /*
   *                        ---- MultiQC ----
   *
   * MultiQC agrega reportes de FastQC, Fastp y otras herramientas en un
   * unico reporte HTML interactivo.
   *
   * Input:
   *   - qc_files (path): Coleccion de reportes (ZIPs de FastQC, JSONs de Fastp, etc.).
   *   - outdir (val): Ruta de salida del reporte.
   *
   * Output:
   *   - multiqc_report.html: Reporte HTML interactivo de MultiQC.
   */
  cache 'lenient'
  container 'ariel-env:latest'
  publishDir { outdir }, mode: 'copy'

  input:
    path qc_files
    val outdir

  output:
    path("multiqc_report.html")

  script:
  """
    multiqc .
  """
}
