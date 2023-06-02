#fasttree for genotube

FROM ubuntu:latest

RUN apt-get update && apt-get install -y  \
    fasttree
