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

# #Install Miniconda
RUN cd $WD && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p $WD/miniconda
RUN $WD/miniconda/bin/conda init bash
RUN $WD/miniconda/bin/conda install -y pip pandas numpy matplotlib seaborn
RUN $WD/miniconda/bin/pip install 'xgboost==2.0.3' 'scikit-learn==1.4.2'
#RUN $WD/miniconda/bin/conda install -y -c conda-forge scikit-learn py-xgboost

CMD ["/bin/bash"]
