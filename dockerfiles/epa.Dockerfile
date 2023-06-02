#epa-ng for genotube

FROM ubuntu:latest as builder

RUN apt-get update && apt-get install -y  \
    wget make g++ cmake flex bison zlib1g-dev

RUN wget https://github.com/Pbdas/epa-ng/archive/v0.3.8.tar.gz

RUN tar -xzvf v0.3.8.tar.gz

RUN cd epa-ng-0.3.8 &&  \
    make

FROM ubuntu:latest as base

RUN apt-get update && apt-get install -y libgomp1

COPY --from=builder epa-ng-0.3.8/bin/epa-ng /bin/epa-ng
