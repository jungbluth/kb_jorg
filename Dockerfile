FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get install -y wget git

RUN wget https://github.com/bachev/mira/releases/download/V5rc1/mira_V5rc1_linux-gnu_x86_64_static.tar.bz2 && \
    tar -xvf mira_* && \
    rm mira_V5rc1_linux-gnu_x86_64_static.tar.bz2

RUN conda install -c bioconda seqtk

RUN conda install -c bioconda bwa

RUN conda install -c bioconda last && touch /var/log/wtmp

RUN mkdir pilon-1.2.3 && \
    wget https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar -P pilon-1.2.3

RUN wget http://eddylab.org/infernal/infernal-1.1.3-linux-intel-gcc.tar.gz && \
    tar -xvf infernal-* && \
    rm infernal-1.1.3-linux-intel-gcc.tar.gz

RUN git clone https://github.com/lmlui/Jorg && \
    chmod +x /Jorg/jorg

RUN cd /Jorg/Example && \
    wget https://zenodo.org/record/3889132/files/SRX3307784_clean.fastq.gz && \
    gunzip SRX3307784_clean.fastq.gz && \
    mv SRX3307784_clean.fastq bin.186.fastq

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

ENV PATH=/mira_V5rc1_linux-gnu_x86_64_static/bin:$PATH
ENV PATH=/pilon-1.2.3:$PATH
ENV PATH=/infernal-1.1.3-linux-intel-gcc/binaries:$PATH
ENV PATH=/Jorg:$PATH

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
