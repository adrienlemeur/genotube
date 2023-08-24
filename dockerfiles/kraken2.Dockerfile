#kraken2 for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y kraken2
