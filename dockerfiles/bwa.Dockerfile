#bwa for genotube v0.0.1

FROM ubuntu:latest as builder

RUN apt-get update && apt-get install -y wget bzip2

RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2

RUN tar jxf bwa-mem2-2.0pre2_x64-linux.tar.bz2

FROM ubuntu:latest as base

RUN apt-get update && apt-get install -y samtools

COPY --from=builder bwa-mem2-2.0pre2_x64-linux/bwa* /bin