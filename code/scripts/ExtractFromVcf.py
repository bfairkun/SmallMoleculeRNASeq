#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : readVcf
# @created     : Wednesday Dec 15, 2021 20:25:22 CST
#
# @description : 
######################################################################

import sys
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "ClinVar/clinvar.pangolin.vcf" ,"scratch/test.txt.gz"]

_, VcfIn, TsvOut = sys.argv

import vcf
import collections
import gzip



vcf_reader = vcf.Reader(open(VcfIn, 'r'))
with gzip.open(TsvOut, 'wt') as f:
    for record in vcf_reader:
        try:
            OmimOut = ' '.join([i for i in record.INFO["CLNDISDB"] if "OMIM" in i])
            _ = f.write(f'{record.CHROM}\t{record.POS}\t{record.INFO["CLNSIG"][0]}\t{record.INFO["CLNDN"]}\t{record.INFO["MC"][0]}\t{record.INFO["Pangolin"][0]}\t{OmimOut}\n')
        except:
            pass

# vcf_reader = vcf.Reader(open("ClinVar/clinvar.pangolin.vcf", 'r'))
# for i, record in enumerate(vcf_reader):
#     # if i < 1000:
#     try:
#         OmimOut = [i for i in record.INFO["CLNDISDB"] if "OMIM" in i]
#         # if len(OmimOut) > 1:
#         print(OmimOut)
#     except: pass
