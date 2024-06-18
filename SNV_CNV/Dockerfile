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

# Updating ubuntu and installing other necessary software
RUN apt-get update --yes \
    && apt-get install --yes gcc python3-dev python3-setuptools ca-certificates curl gnupg lsb-release wget build-essential make libssl-dev zlib1g-dev git-all python3-pip nodejs unzip rsync r-base r-base-dev \
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

RUN apt-get install --yes apt-transport-https ca-certificates gnupg curl \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 871920D1991BC93C --keyring /usr/share/keyrings/cloud.google.gpg add -

RUN apt-get update && apt-get install google-cloud-cli google-cloud-cli-anthoscli

WORKDIR /opt

RUN apt-get install -y openjdk-17-jdk openjdk-17-jre

RUN git clone https://github.com/broadinstitute/picard.git

WORKDIR /opt/picard

RUN ./gradlew shadowJar

WORKDIR /opt

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip \
    && unzip gatk-4.1.4.1.zip
WORKDIR /opt/gatk-4.1.4.1

# WORKDIR /opt
# Install Miniconda on a MacM1:
# RUN cd $WD && wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-aarch64.sh && bash Miniconda3-py39_4.9.2-Linux-aarch64.sh -b -p $WD/miniconda

# In a Linux or MacIntel:
RUN cd $WD && wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh && bash Miniconda3-4.7.12-Linux-x86_64.sh -b -p $WD/miniconda
RUN $WD/miniconda/bin/conda init bash
RUN $WD/miniconda/bin/conda update -n base -c defaults conda
RUN $WD/miniconda/bin/conda env export -n gatk -f gatkcondaenv.yml

CMD ["/bin/bash"]
RUN conda activate gatk