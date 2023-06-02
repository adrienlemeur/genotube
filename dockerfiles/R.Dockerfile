#R for genotube v0.0.1

FROM ubuntu:latest

RUN apt-get update 	&& apt-get install -y --no-install-recommends \
    ca-certificates ed fonts-texgyre less locales  \
    vim-tiny wget sudo && \
    rm -rf /var/lib/apt/lists/*

ENV CONTAINER_TIMEZONE Europe/Paris
ENV TZ Europe/Paris

RUN sudo echo "Europe/Paris" > /etc/timezone
RUN echo "fr_FR.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen fr_FR.UTF-8 \
  && /usr/sbin/update-locale LANG=fr_FR.UTF-8 \

ENV LC_ALL fr_FR.UTF-8
ENV LANG fr_FR.UTF-8

RUN apt-get update &&  \
    apt-get install -y --no-install-recommends  \
    gcc-9-base libopenblas0-pthread littler r-cran-littler  \
    r-base