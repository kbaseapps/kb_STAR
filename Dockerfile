FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install coverage

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

###### STAR installation                                                                                                                                       
#  Directions from https://github.com/alexdobin/STAR 
#  with more details at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/ 
#  Download tarball from https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz, 
#  untar and build STAT.

WORKDIR /kb/module
RUN \
  wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
  tar -zxf 2.5.3a.tar.gz && \
  ln -s STAR-2.5.3a STAR && \
  rm -rf 2.5.3a.tar.gz && \
  cd STAR/source && \
  make && \
  cp STAR /kb/deployment/bin/.  

# The genome directory where the genome indexes are stored. 
# This directory has to be created (with mkdir) before STAR run
# and needs to have writing permissions. 
# The file system needs to have at least 100GB of disk space available for a typical mammalian genome. 
RUN mkdir -p /kb/module/STAR_genome_dir
RUN chmod -R a+rw /kb/module/STAR_genome_dir

#WORKDIR /kb/module
#RUN \
  #mkdir /STAR_genome_directory && \
  #chmod -R a+rw /STAR_genome_directory

# -----------------------------------------
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
