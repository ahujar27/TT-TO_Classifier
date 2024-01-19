FROM ubuntu

#When running this image, mount the docker socket so the workflow can launch docker containers from within this container
#docker run -it -v /var/run/docker.sock:/var/run/docker.sock ubuntu

LABEL maintainer="m.inkman@wustl.edu"
LABEL description="A docker container to run the TUO classifier and prep input data for the classifier using the GDC DNA-seq mutation and copy number pipelines"

#Set environmental variables
ENV WD=/opt \
DEBIAN_FRONTEND=noninteractive

WORKDIR $WD
ADD . $WD

#Updating ubuntu and installing other necessary software
RUN apt-get update --yes \
&& apt-get install --yes ca-certificates curl gnupg lsb-release wget build-essential make libssl-dev zlib1g-dev git python3-pip nodejs unzip rsync r-base r-base-dev \
&& apt-get clean all

# Add source to install docker via apt-get
RUN install -m 0755 -d /etc/apt/keyrings \
 && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg \
 && chmod a+r /etc/apt/keyrings/docker.gpg
    
RUN echo \
 "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
 "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
 tee /etc/apt/sources.list.d/docker.list > /dev/null

# Install docker
RUN apt-get update --yes \
&& apt-get install --yes docker-ce-cli

RUN mkdir references

WORKDIR /opt/references/

RUN wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O GRCh38.d1.vd1.fa.tar.gz \
&& tar xf GRCh38.d1.vd1.fa.tar.gz \
&& rm GRCh38.d1.vd1.fa.tar.gz \
&& wget https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225 \
&& tar xf 25217ec9-af07-4a17-8db9-101271ee7225 \
&& rm 25217ec9-af07-4a17-8db9-101271ee7225 \
&& wget https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5 \
&& tar xf 2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5 \
&& rm 2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5

# RUN apt-get install --yes apt-transport-https ca-certificates gnupg curl \
#     && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
#     && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
    
# RUN apt-get update && apt-get install google-cloud-cli

RUN apt-get install -y openjdk-17-jdk openjdk-17-jre

#Install picard tools

RUN git clone https://github.com/broadinstitute/picard.git

WORKDIR picard
RUN ./gradlew shadowJar

WORKDIR $WD 

RUN apt-get -y install libbz2-dev liblzma-dev libncurses5-dev tabix 

RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
&& tar xjvf samtools-1.18.tar.bz2 \
&& cd samtools-1.18 \
&& ./configure \
&& make \
&& make install

#Install Miniconda
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip \
&& unzip gatk-4.4.0.0.zip
WORKDIR gatk-4.4.0.0
# #Install Miniconda
RUN cd $WD && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p $WD/miniconda
RUN $WD/miniconda/bin/conda init bash
RUN $WD/miniconda/bin/conda env create -f gatkcondaenv.yml
RUN $WD/miniconda/bin/conda install -n gatk -qy -c conda-forge -c bioconda -c defaults ensembl-vep==102.0 htslib==1.10.2 bcftools==1.10.2 samtools==1.10 ucsc-liftover==377

WORKDIR $WD

RUN export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4` \
&& curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-* \
&& perl vcf2maf.pl --man \
&& perl maf2maf.pl --man

RUN sed -i 's/\$ENV{HOME}/$WD/g' /opt/mskcc-vcf2maf-754d68a/vcf2maf.pl
RUN sed -i '15i\my ($WD) = ("/opt");' /opt/mskcc-vcf2maf-754d68a/vcf2maf.pl

RUN mkdir -p $WD/.vep/homo_sapiens/102_GRCh38/ \
&& rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh38.tar.gz $WD/.vep/ \
&& tar -zxf $WD/.vep/homo_sapiens_vep_102_GRCh38.tar.gz -C $WD/.vep/ \
&& rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/fasta/homo_sapiens/dna_index/ $WD/.vep/homo_sapiens/102_GRCh38/

WORKDIR $WD


CMD ["/bin/bash"]
