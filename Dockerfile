FROM ubuntu:16.04
LABEL maintainer "Rico Derks" r.j.e.derks@lumc.nl

# install stuff for R
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    wget \
  && . /etc/environment

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran40/"

RUN gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add -

# For newer versions
#RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# install R
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libxml2-dev \
    r-base \
    r-base-dev \
  && . /etc/environment

# install some R packages
RUN R -e 'install.packages("XML")'

# install other stuff
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gnuplot \
    apache2 \
    apache2-utils \
    php \
    libapache2-mod-php \
    libgd-dev \
  && . /etc/environment

WORKDIR /var/www/html

# copy all stuff for msalign2
ADD ./msalign2.tar.gz /var/www/html

RUN chmod +x ./*.sh

# compile msalign2
RUN gcc -o msalign2 base64.c ramp.c msalign2.c -I. -lgd -lm -lz -std=gnu99

# some folders
RUN mkdir -p /var/www/html/data
RUN mkdir -p /var/www/html/temp

EXPOSE 80

CMD ["apache2ctl", "-D","FOREGROUND"]
