FROM kbase/sdkbase2:latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN pip install pathos

###### STAR installation                                                                                                                                       
#  Directions from https://github.com/alexdobin/STAR 
#  with more details at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/ 
#  Download tarball from https://github.com/alexdobin/STAR/archive/2.6.1a.tar.gz, 
#  untar and build STAR.
ENV VERSION=2.6.1a

WORKDIR /kb/module
RUN wget https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz && \
  tar -zxf ${VERSION}.tar.gz && \
  ln -s STAR-${VERSION} STAR && \
  rm -rf ${VERSION}.tar.gz && \
  cd STAR/source && \
  make && \
  cp STAR /kb/deployment/bin/.  

# The genome directory where the genome indexes are stored. 
# This directory has to be created (with mkdir) before STAR run
# and needs to have writing permissions. 
# The file system needs to have at least 100GB of disk space available for a typical mammalian genome. 
# RUN mkdir -p /kb/module/STAR_Genome_index
# RUN chmod -R a+rw /kb/module/STAR_Genome_index

# -----------------------------------------
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
