#fastq_screen for genotube v0.0.1

FROM ubuntu:latest as builder

RUN apt-get update && \
    apt-get install -y wget tar

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.14.0.tar.gz

RUN tar xvzf fastq_screen_v0.14.0.tar.gz

FROM ubuntu:latest as base

RUN apt-get update && \
    apt-get install -y perl

COPY --from=builder fastq_screen_v0.14.0/fastq_screen /bin/fastq_screen

COPY --from=builder fastq_screen_v0.14.0/interactive_graphs.js /bin/interactive_graphs.js