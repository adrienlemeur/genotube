#bbtools for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    wget tar openjdk-8-jre-headless

RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz

RUN tar -xvzf BBMap_39.01.tar.gz

RUN mv bbmap/* /bin/