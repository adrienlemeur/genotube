#python for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update && \
     apt-get install -y python2

RUN mv bin/python2 bin/python