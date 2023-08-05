#!/usr/bin/env bash

binID=${1}

seqtk seq -1 -C ${binID}.fastq | seqtk rename - ${binID}_ | sed "s/\(^@${binID}_[0-9][0-9]*$\)/\1\/1/" > ${binID}_1.fastq
seqtk seq -2 -C ${binID}.fastq | seqtk rename - ${binID}_ | sed "s/\(^@${binID}_[0-9][0-9]*$\)/\1\/2/" > ${binID}_2.fastq