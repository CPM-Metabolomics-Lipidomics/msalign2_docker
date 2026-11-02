FROM ubuntu:16.04
LABEL maintainer="Rico Derks <r.j.e.derks@lumc.nl>"

# install R and the XML package from the Ubuntu repository
# (the CRAN xenial-cran40 repository no longer installs on 16.04: r-base-dev needs libdeflate-dev)
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    r-base \
    r-cran-xml \
  && . /etc/environment

# install other stuff
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gnuplot \
    apache2 \
    apache2-utils \
    php \
    libapache2-mod-php \
    libgd-dev \
    gcc \
    libc6-dev \
  && . /etc/environment

WORKDIR /var/www/html

# copy all stuff for msalign2
ADD ./msalign2.tar.gz /var/www/html

RUN chmod +x ./*.sh

# compile msalign2
RUN gcc -o msalign2 base64.c ramp.c msalign2.c -I. -lgd -lm -lz -std=gnu99

# add some folders
RUN mkdir -p /var/www/html/data
RUN mkdir -p /var/www/html/temp

# make writable for the webserver
RUN chmod a+rwx /var/www/html/data
RUN chmod a+rwx /var/www/html/temp

EXPOSE 80

CMD ["apache2ctl", "-D","FOREGROUND"]
