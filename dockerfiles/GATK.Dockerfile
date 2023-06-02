#gatk for genotube v0.0.1

FROM ubuntu:latest as builder

RUN apt-get update && apt-get install -y \
    wget unzip python2 openjdk-8-jdk

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip

RUN unzip gatk-4.3.0.0.zip

FROM ubuntu:latest as base

RUN apt-get update && apt-get install -y \
    python2 openjdk-8-jdk

RUN mv bin/python2 bin/python

COPY --from=builder gatk-4.3.0.0/gatk /bin/gatk

COPY --from=builder gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar /bin/gatk-package-4.3.0.0-local.jar
