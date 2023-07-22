FROM kbase/sdkpython:3.8.10
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get install -y wget git

RUN apt-get update && apt-get install -y build-essential make flex libexpat1-dev libboost-all-dev xxd zlib1g-dev bowtie2 bedtools

RUN wget https://github.com/bachev/mira/releases/download/V5rc2/mira-V5rc2.tar.bz2 && \
    tar -xvf mira-* && \
    cd mira-V5rc2 && \
    ./configure && \
    make && \
    make install && \
    rm /mira-V5rc2.tar.bz2

RUN git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    make && \
    make install

RUN conda install -c bioconda bwa

RUN conda install -c bioconda last && touch /var/log/wtmp

RUN conda install -c bioconda seqkit

RUN conda install -c bioconda circos

RUN conda install pandas

#RUN mkdir pilon-1.2.3 && \
#    wget https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar -P pilon-1.2.3

# RUN wget http://eddylab.org/infernal/infernal-1.1.3-linux-intel-gcc.tar.gz && \
#    tar -xvf infernal-* && \
#    rm infernal-1.1.3-linux-intel-gcc.tar.gz

RUN git clone https://github.com/jungbluth/Jorg && \
    chmod +x /Jorg/jorg

# RUN mv SRX3307784_clean.fastq.gz /Jorg/Example

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

ENV PATH=/kb/module/lib/kb_jorg/bin:$PATH
ENV PATH=/kb/module/lib/kb_jorg/bin/bbmap:$PATH
ENV PATH=/kb/module/lib/kb_jorg/bin/samtools/bin:$PATH
ENV PATH=/pilon-1.2.3:$PATH
ENV PATH=/infernal-1.1.3-linux-intel-gcc/binaries:$PATH
ENV PATH=/Jorg:$PATH

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
