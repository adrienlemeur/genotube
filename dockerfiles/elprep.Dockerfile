#elprep for genotube v0.0.1

FROM ubuntu:latest as builder

RUN apt-get update && \
    apt-get install -y wget tar

RUN wget https://github.com/ExaScience/elprep/releases/download/v5.1.2/elprep-v5.1.2.tar.gz

RUN tar xvf elprep-v5.1.2.tar.gz

FROM ubuntu:latest as base

COPY --from=builder elprep /bin/