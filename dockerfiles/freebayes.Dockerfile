#freebayes for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y freebayes tabix