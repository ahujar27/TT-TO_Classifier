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
    && apt-get install --yes gcc python3-dev python3-setuptools ca-certificates curl gnupg lsb-release wget build-essential make libssl-dev zlib1g-dev git python3-pip nodejs unzip rsync r-base r-base-dev \
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

# RUN mkdir references

# WORKDIR /opt/references/

# RUN wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O GRCh38.d1.vd1.fa.tar.gz \
# && tar xf GRCh38.d1.vd1.fa.tar.gz \
# && rm GRCh38.d1.vd1.fa.tar.gz \
# && wget https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225 \
# && tar xf 25217ec9-af07-4a17-8db9-101271ee7225 \
# && rm 25217ec9-af07-4a17-8db9-101271ee7225 \
# && wget https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5 \
# && tar xf 2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5 \
# && rm 2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5

# RUN apt-get install --yes apt-transport-https ca-certificates gnupg curl \
#     && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
#     && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

# RUN apt-get update && apt-get install google-cloud-cli

# RUN gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf . \
# && gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx .

# RUN gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz . \
# && gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi .

# RUN gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz . \
# && gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi .

# RUN gsutil cp gs://gatk-test-data/mutect2/Homo_sapiens_assembly38.index_bundle .

# WORKDIR /opt/

# RUN apt-get update \
# && apt-get install -y autoconf automake libtool pkg-config g++ make libboost-dev

# # RUN git clone https://gitlab.com/german.tischler/libmaus2

# # WORKDIR /opt/libmaus2

# # RUN libtoolize \
# # && aclocal \
# # && autoreconf -i -f \
# # && ./configure \
# # && make \
# # && make install

# # WORKDIR /opt

# # RUN git clone https://gitlab.com/german.tischler/biobambam2.git

# # WORKDIR /opt/biobambam2

# # RUN autoreconf -i -f \
# # && ./configure --with-libmaus2=/usr/local --prefix=/opt/biobambam2 \
# # && make \
# # && make install \
# # && ldconfig

# # WORKDIR /opt

# RUN apt-get -y install libbz2-dev liblzma-dev libncurses5-dev tabix

# # RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
# # && tar xjvf samtools-1.18.tar.bz2 \
# # && cd samtools-1.18 \
# # && ./configure \
# # && make \
# # && make install

# # RUN git clone https://github.com/lh3/bwa.git

# # WORKDIR /opt/bwa

# # RUN make

# # # RUN /opt/bwa/bwa index /opt/references/GRCh38.d1.vd1.fa

WORKDIR /opt

RUN apt-get install -y openjdk-17-jdk openjdk-17-jre

RUN git clone https://github.com/broadinstitute/picard.git

WORKDIR /opt/picard

RUN ./gradlew shadowJar

WORKDIR /opt

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip \
    && unzip gatk-4.4.0.0.zip
WORKDIR /opt/gatk-4.4.0.0

# WORKDIR /opt
# #Install Miniconda
RUN cd $WD && wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-aarch64.sh && bash Miniconda3-py39_4.9.2-Linux-aarch64.sh -b -p $WD/miniconda
RUN $WD/miniconda/bin/conda init bash
RUN $WD/miniconda/bin/conda update -n base -c defaults conda
RUN $WD/miniconda/bin/conda env export -n gatk -f gatkcondaenv.yml

CMD ["/bin/bash"]
