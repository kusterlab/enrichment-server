#How to start shell inside docker: docker run -it --entrypoint /bin/bash enrichment_server
#How to run tests: docker run -it --entrypoint /bin/sh enrichment_server -c "Xvfb :1 -screen 0 1024x768x24 & poetry run pytest -v"
#How to run it as a server: docker run --network host enrichment_server

#We need R 4.3.2
FROM rocker/r-ver:4.3.2

MAINTAINER Julian Müller "julian2.mueller@tum.de"

WORKDIR /app

# Tell docker that we don't want to be bothered with questions
ARG DEBIAN_FRONTEND=noninteractive


RUN apt-get update && apt-get install -y \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libomp5 \
    libopenblas-dev \
    libfontconfig1-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    openjdk-17-jre \
    xvfb \
    git \
    python3.11 \
    && rm -rf /var/lib/apt/lists/*
#Get Poetry
RUN curl -sSL 'https://install.python-poetry.org' | python3.11 -
#Add it to Path for the docker
ENV PATH="${PATH}:/root/.local/bin:/root/.poetry/bin"
#Add it to Path in case someone logs into an interactive shell
RUN echo 'export PATH="/root/.local/bin:$PATH"' >> ~/.bashrc


#Ubuntu 22.04 only support libssl3. as we need libssl1.1 we need to install it separately
#Note that we manually need to update this url when this package is updated in the ubuntu archives, e.g. ubuntu2.22 => ubuntu2.23
RUN curl http://nz2.archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2.23_amd64.deb --output libssl.deb
RUN dpkg -i libssl.deb
RUN rm libssl.deb

#Set Up PTM-SEA
RUN git clone -b master --single-branch https://github.com/broadinstitute/ssGSEA2.0.git
WORKDIR /app/ssGSEA2.0
RUN git checkout 4b5198fb5a19759eec2625ba5dde03fc861d96ac
#Fix a sloppy line from ssGSEA by hard-coding the script directory
RUN sed -i '8,13d; 6i script.dir<-"../ssGSEA2.0"' ssgsea-cli.R

#Set Up RoKAI
WORKDIR /app
RUN git clone -b master --single-branch https://github.com/serhan-yilmaz/RokaiApp.git
WORKDIR /app/RokaiApp
RUN git checkout e2eb10a72c6dc5e19572e1013313fb4b7a1310ab

# Copy the current directory contents into the container at /app
COPY db /app/db
COPY fixtures /app/fixtures
COPY CPLEX /app/CPLEX
COPY Cytoscape_v3.10.1 /app/Cytoscape_v3.10.1
COPY CytoscapeConfiguration /root/CytoscapeConfiguration

#Install R-related stuff
COPY flask_server/renv /app/flask_server/renv
COPY flask_server/.renvignore /app/flask_server/.renvignore
COPY flask_server/.Rprofile /app/flask_server/.Rprofile
COPY flask_server/renv.lock /app/flask_server/renv.lock
WORKDIR /app/flask_server
RUN Rscript -e "install.packages('renv')"
RUN Rscript -e "renv::restore()"
#Stringi is broken in the renv for some reason, reinstall it (else it is missing libicui18n.so.66)
RUN Rscript -e "install.packages('stringi')"

#Install Python-related stuff
COPY flask_server/poetry.lock /app/flask_server/poetry.lock
COPY flask_server/pyproject.toml /app/flask_server/pyproject.toml
WORKDIR /app
RUN poetry install -C flask_server

WORKDIR /app/flask_server
#The KSTAR Python Package needs some initial setting up
RUN poetry run python -c "from kstar import config; \
    config.install_resource_files(); \
    config.NETWORK_DIR, config.NETWORK_Y_PICKLES, config.NETWORK_ST_PICKLES = config.update_network_directory('../db/kstar/NetworKIN')"

#Copy rest of the flask_server directory including the modules and tests
COPY flask_server /app/flask_server

#Tell Cytoscape to use the Xvfb virtual display
ENV DISPLAY=:1

#Start Xvfb for headless Cytoscape, then run the python server
ENTRYPOINT Xvfb :1 -screen 0 1024x768x24 & poetry run python enrichment_server.py

EXPOSE 4321
