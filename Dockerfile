#How to start shell inside docker: docker run -it --entrypoint /bin/bash enrichment_server
#How to run tests: docker run -it enrichment_server poetry run pytest -v
#How to run it as a server: docker run --network host enrichment_server #TODO: It feels like I should not have to use '--network host', but I get a connection refused error if I don't.

#TODO: Find a slimmer image
FROM ubuntu:24.04

MAINTAINER Julian Müller "julian2.mueller@tum.de"

WORKDIR /app

# Tell docker that we don't want to be bothered with questions
ARG DEBIAN_FRONTEND=noninteractive


RUN apt-get update
RUN apt-get install -y curl libcurl4-openssl-dev libssl-dev libomp5 libopenblas-dev libfontconfig1-dev \
    libxml2-dev libharfbuzz-dev libfribidi-dev libtiff5-dev libssl-dev python3-poetry openjdk-17-jre xvfb git r-base #Will install R=4.3.2
#Make sure python3.11 is installed, the gitlab pipeline would sometimes default to python3.12
RUN apt-get install -y python3.11

##Set Up PTM-SEA
RUN git clone -b master --single-branch https://github.com/broadinstitute/ssGSEA2.0.git
WORKDIR /app/ssGSEA2.0
RUN git checkout 4b5198fb5a19759eec2625ba5dde03fc861d96ac

#Set Up RoKAI
WORKDIR /app
RUN git clone -b master --single-branch https://github.com/serhan-yilmaz/RokaiApp.git
WORKDIR /app/RokaiApp
RUN git checkout e2eb10a72c6dc5e19572e1013313fb4b7a1310ab

# Copy the current directory contents into the container at /app
COPY flask_server /app/flask_server
COPY db /app/db
COPY fixtures /app/fixtures
COPY CPLEX /app/CPLEX
COPY Cytoscape_v3.10.1 /app/Cytoscape_v3.10.1
COPY CytoscapeConfiguration /root/CytoscapeConfiguration

#Install R-related stuff
RUN Rscript -e "install.packages('packrat')"
RUN Rscript -e "packrat::restore('/app/flask_server')"

WORKDIR /app
#Install Python-related stuff
RUN poetry install -C flask_server


WORKDIR /app/flask_server

CMD ["poetry", "run", "python", "enrichment_server.py"]

EXPOSE 4321
