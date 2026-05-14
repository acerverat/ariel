#!/bin/bash

set -e

if [[ -z "$1" ]]; then
	echo "No se introdujo una ruta"
	echo "Para generar el directorio de referencias es necesario correr el script de la siguiente manera:"
	echo " bash generaReferencias.sh /la_ruta_completa"
	echo "Ejemplo: bash generaReferencias.sh /home/yun/Documentos/referencias"
	exit 1
fi

rutaReferencias="$1"
echo "la ruta para almacenar las referencias es: $rutaReferencias"

mkdir -p "$rutaReferencias"
cd "$rutaReferencias"
mkdir -p GRCh38_no_alt
cd GRCh38_no_alt


# descarga referencias del repositorio de Github https://github.com/stjude/CICERO#reference
wget -c https://zenodo.org/records/5088371/files/reference.tar.gz

tar -xvzf reference.tar.gz

#modificamos la estructura de las referencias, eliminando el nivel reference
mv reference cicero_references/
rm reference.tar.gz

# descarga anotaciones
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gff3.gz

# descomprime
gunzip -f gencode.v42.annotation.gtf.gz
gunzip -f gencode.v42.annotation.gff3.gz




# rutas
fasta=$PWD/cicero_references/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa
gtf=$PWD/gencode.v42.annotation.gtf
gff3=$PWD/gencode.v42.annotation.gff3

# Indice STAR

mkdir -p STAR_2.7.10b_index

docker run  -u $(id -u):$(id -g) --rm \
		-v $PWD:$PWD \
	        -w $PWD \
		ariel-env \
		STAR \
		--runMode genomeGenerate \
		--genomeDir STAR_2.7.10b_index \
		--genomeFastaFiles $fasta \
		--sjdbGTFfile $gtf \
		--runThreadN 4

# Base de datos de FusionCatcher
# Se descarga usando el contenedor para garantizar compatibilidad con la versión instalada.
# El directorio fusioncatcher_db se monta en /opt/fusioncatcher/data/current dentro del
# contenedor cada vez que se ejecuta el proceso FusionCatcher.

mkdir -p fusioncatcher_db

docker run -u $(id -u):$(id -g) --rm \
    -v $PWD/fusioncatcher_db:/opt/fusioncatcher/data/current \
    ariel-env \
    bash -c "/opt/fusioncatcher/bin/download-human-db.sh"

# Índice de Salmon (GRCh38, partial SA index) desde refgenie
wget -c "http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_partial_sa_index?tag=default" \
    -O salmon_index.tar.gz

tar -xzf salmon_index.tar.gz

# refgenie extrae a salmon_partial_sa_index/default/
mv salmon_partial_sa_index/default salmon_index
rm -rf salmon_partial_sa_index salmon_index.tar.gz
