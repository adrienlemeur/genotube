#raxml for genotube

FROM ubuntu:latest as builder

RUN apt-get update && apt-get install -y  \
    git flex bison libgmp3-dev cmake wget unzip

RUN wget https://github.com/amkozlov/raxml-ng/releases/download/0.9.0/raxml-ng_v0.9.0_linux_x86_64.zip && \
  unzip raxml-ng_v0.9.0_linux_x86_64.zip

FROM ubuntu:latest as base

RUN mkdir raxml_ng

COPY --from=builder raxml-ng /raxml_ng/raxml-ng

ENV PATH="${PATH}:/standard-RAxML-8.2.12:/raxml_ng"
