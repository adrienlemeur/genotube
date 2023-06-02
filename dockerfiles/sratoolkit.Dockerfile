#sra-toolkit for genotube v0.0.1

FROM ubuntu:latest as builder

RUN apt-get update && apt-get install -y \
    wget pigz

RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz

RUN tar xzvf sratoolkit.3.0.0-ubuntu64.tar.gz

FROM ubuntu:latest as base

RUN apt-get update && apt-get install -y \
    wget pigz

COPY --from=builder sratoolkit.3.0.0-ubuntu64/bin /bin
