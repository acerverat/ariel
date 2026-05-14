FROM ubuntu:24.04

# Build the image: docker build -t ariel-env .


ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    ca-certificates \
    cmake \
    curl \
    default-jdk \
    gawk \
    gcc \
    git \
    gmap \
    gzip \
    g++ \
    libcurl4-openssl-dev \
    libssl-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    r-base \
    r-base-dev \
    r-cran-devtools \
    tar \
    unzip \
    wget \
    zip \
    zlib1g \
    zlib1g-dev \
  && apt-get clean

RUN apt-get install -y r-bioc-genomicalignments r-bioc-genomicranges 
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository --update -y ppa:ubuntu-toolchain-r/test 
RUN apt-get update -y 
RUN apt-get -y --fix-broken install gcc-11 g++-11


RUN Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'install.packages("ranger")'
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("ranger")' 
RUN Rscript -e 'install.packages("tidyverse",dependencies=T)'
RUN Rscript -e 'BiocManager::install("tximport")'
RUN Rscript -e 'install.packages("readr")'

WORKDIR /usr/local/lib

###### Samtools & friends #####
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
RUN git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    make
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar xjvf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    ./configure --prefix=/usr/local/ && \
    make && \
    make install

####### STAR-Fusion #######
RUN wget -qO - 'https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz' | \
tar --strip-components=3 -C /usr/local/bin -xzf - 'STAR-2.7.10b/bin/Linux_x86_64/STAR'


##### Salmon #####

RUN curl -L https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz | tar xzv && \
    ln -s /usr/local/lib/salmon-1.9.0_linux_x86_64/bin/salmon /usr/local/bin/salmon


# install arriba
RUN wget -qO - 'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' | tar -xzf - --exclude='arriba*/.git'

# make wrapper script for download_references.sh
RUN echo '#!/bin/bash\n\
cd /references\n\
/arriba*/download_references.sh $1 && \\\n\
ASSEMBLY=$(sed -e "s/viral+.*//" -e "s/+.*//" <<<"$1") && \\\n\
cp /arriba*/database/*$ASSEMBLY* /references' > /usr/local/bin/download_references.sh && \
chmod a+x /usr/local/bin/download_references.sh

# make wrapper script for run_arriba.sh
RUN echo '#!/bin/bash\n\
cd /output && \\\n\
/arriba*/run_arriba.sh /references/STAR_index_* /references/*.gtf /references/*.fa /references/blacklist_*.tsv.gz /references/known_fusions_*.tsv.gz /references/protein_domains_*.gff3 ${THREADS-8} /read1.fastq.gz $(ls /read2.fastq.gz 2> /dev/null)' > /usr/local/bin/arriba.sh && \
chmod a+x /usr/local/bin/arriba.sh

# make wrapper script for draw_fusions.R
RUN echo '#!/bin/bash\n\
Rscript /arriba*/draw_fusions.R --annotation=$(ls /references/*.gtf) --fusions=/fusions.tsv --output=/output/fusions.pdf --proteinDomains=$(ls /references/protein_domains_*.gff3) --alignments=/Aligned.sortedByCoord.out.bam --cytobands=$(ls /references/cytobands_*.tsv)' > /usr/local/bin/draw_fusions.sh && \
chmod a+x /usr/local/bin/draw_fusions.sh

# FusionCatcher
RUN apt-get install -y --no-install-recommends make automake pigz parallel libncurses5-dev libc6-dev zlib1g zlib1g-dev libtbb12 libtbb-dev libtbbmalloc2

RUN wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz && \
    tar -xzf Python-2.7.18.tgz && \
    cd Python-2.7.18 && \
    ./configure --prefix=/usr/local/python2 --enable-optimizations && \
    make -j"$(nproc)" && \
    make install

# Crear symlinks globales
RUN ln -sf /usr/local/python2/bin/python2 /usr/bin/python2 && \
    ln -sf /usr/local/python2/bin/python2 /usr/bin/python && \
    true

# Instalar pip2
RUN curl -s https://bootstrap.pypa.io/pip/2.7/get-pip.py | python2

# Instalar numpy para Python 2
RUN python2 -m pip install numpy
RUN python2 -m pip install xlrd openpyxl
RUN python2 -m pip install biopython==1.70

# Instalar FusionCatcher desde bootstrap (instala en /opt/fusioncatcher)
RUN wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py \
    && python2 bootstrap.py -t --download -y

# Crear las carpetas necesarias y hacer symlink a los datos
RUN mkdir -p /usr/local/lib/fusioncatcher/data && \
    ln -s /opt/fusioncatcher/data/current /usr/local/lib/fusioncatcher/data/current

# Crear archivo de configuración
RUN echo 'FUSIONCATCHER_HOME=/opt/fusioncatcher \nFUSIONCATCHER_DB=/usr/local/lib/fusioncatcher/data/current' > myConfig.txt


# Instalar Fungi
RUN git clone https://alejandra_cervera@bitbucket.org/alejandra_cervera/fungi.git
RUN cd fungi && git fetch && git checkout v2 
RUN cd fungi/resources && ./download.sh 
ARG Fungi=unknown
RUN apt install python3-pandas -y
RUN cd /usr/local/lib/fungi && git pull
RUN ln -s /usr/local/lib/fungi/bin/* /usr/local/bin/
RUN cd /opt/fusioncatcher/bin/
RUN wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/refs/heads/master/bin/label_fusion_genes_znf.py
RUN cp label_fusion_genes_znf.py /opt/fusioncatcher/bin/
RUN chmod +x /opt/fusioncatcher/bin/label_fusion_genes_znf.py

ENV TZ="America/Mexico_City"

CMD ["bash", "-l"]
