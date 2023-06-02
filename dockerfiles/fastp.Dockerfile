#fastp for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y \
    wget

RUN mkdir fastp &&  \
    cd fastp &&  \
    wget http://opengene.org/fastp/fastp.0.23.2 &&  \
    mv fastp.0.23.2 fastp &&  \
    chmod a+x ./fastp &&  \
    mkdir /data # buildkit

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/fastp/

WORKDIR /data