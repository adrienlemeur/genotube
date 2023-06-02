#bracken for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y git make g++ python2

RUN git clone https://github.com/jenniferlu717/Bracken.git

RUN cd Bracken && \
    bash install_bracken.sh && \
    mv bracken /bin/bracken && \
    mkdir /bin/src/ && \
    mv src/* /bin/src/

RUN mv bin/python2 bin/python