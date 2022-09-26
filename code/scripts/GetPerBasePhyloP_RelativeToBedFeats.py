#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : GetPerBasePhyloP_RelativeToBedFeats
# @created     : Monday Sep 26, 2022 09:55:39 CDT
#
# @description : to get per base phyloP scores for each base in input bed file. Each bed feature should be of identicial length.
######################################################################

import sys
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "Arg1" ,"Arg2"]

_, MyArg1, MyArg2 = sys.argv

def main():
    pass

if __name__ == '__main__':
    main()

