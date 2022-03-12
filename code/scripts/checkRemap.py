"""
author: Chao Dai
date: 11/18/21
Description: This script checks the first pass called A-I (A-G, or T-C) editing sites based on 2 critia. It requires 3 input
files: 
    1. a BED file of first pass putative editing sites.
    2. a BAM file from the first alignment. Likely the bam input file for BCFtools mpileup
    3. a BAM file from the second alignment that used first passs putative editing sites as WASP variant file. It should have been
       filtered to include only alignments that overlap with first pass editing sites (vW = 1)

Qualified editing sites must satisfy two checks, implemented in functions checkMapping and checkSite:
    1. reads that aligned to the editing site during the first mapping procedure and the second mapping procedure must be consistent
    2. after realigning the reads using first pass editing sites as variant file, there must be at least 1 A>G or T>C mismatch

"""




def getComplementary(letter):
    """
    Description: return complementary base, case sensitive
    """
    COMPLEMENTARY = {'A':'T', 'C':'G', 'G':'C','T':'A'}
    if letter in COMPLEMENTARY.keys():
        v = COMPLEMENTARY.get(letter)
    else:
        v = 'N'
    return v


def checkMapping(ch, st, bam1, bam2, min_agreement = 0.499):
    """
    Description: compare the 1st and 2nd alignment at given coordinate
    ch: chr
    st: 0-based pos
    bam1: first alignment as pysam alignment file object
    bam2: second alignment (filtered) as pysam alignment file object
    min_agreement: minimum agreement set at 0.499 by default
    
    Returns a tuple check = True/False means passing check; agreement = % of second alignment that agrees with the first alignment
    """
    en = st + 1 # end pos
    aln1 = [] # read names of alignments overlapping with given coordinates from bam1
    for r in bam1.fetch(ch, st, en):
        aln1.append(r.qname)

    aln2 = [] # read names of alignments overlapping with given coordinates from bam2
    for r in bam2.fetch(ch, st, en):
        aln2.append(r.qname)
    if len(aln1) > 0 and len(aln2) > 0:
        # % aligned reads from bam2 that are aligned at the same locus in bam1
        agreement = sum([l in aln1 for l in aln2])/len(aln2)
    else:
        agreement = 0
    if agreement > min_agreement:
        # if agreement > threshold, return that reads are generally mapped at the same locus in both mappings
        check = True
    else:
        check = False
    return check, agreement


def checkSite(chrom, start, bamFile, min_edit = 1):
    """
    chrom: chr
    start: start 0 based
    bamFile: pysam alignment file object
    min_edit: minimum number of reads supporting AG or TC editing
    
    Return a tuple passCheck = True / False, N_A2I = number of supporting mismatch reads
    """
    end = start + 1
    ref_query_base_pair = [] # store ref_query_base_pair
    passCheck = False
    if len(list(bamFile.fetch(chrom, start , end))) > 0: # must have alignments
        for r in bamFile.fetch(chrom, start, end):
            if r.has_tag('MD'): # make sure MD tag exisits
                # pysam.get_aligned_pairs return list of aligned query and ref postions
                # here tup[1] == start selects only the base pair that at given coordinate
                aligned_pair = [ tup for tup in 
                                 r.get_aligned_pairs(matches_only=True, with_seq=True) 
                                 if tup[1] == start
                                ]
                #print(aligned_pair, len(aligned_pair)) # diagnostic
                if len(aligned_pair) > 0: # make sure get_aligned_pairs() return at least 1 result
                    query_base_pos, ref_base_pos, ref_base = aligned_pair[0]
                    ref_base = ref_base.upper()
                    query_base = r.query_sequence[query_base_pos]
                    ref_query_base_pair.append("".join((ref_base, query_base)))
                    #print(r.qname, query_base_pos, ref_base_pos, ref_query_base_pair, passCheck) # diagnostic
                else:
                    #print("no basepair match found")#diagnostic
                    pass

            else: # what happens when there is no MD tag?
                print('No MD') # diagnostic
                pass

        if len(ref_query_base_pair) > 0 and \
                sum([pair in ['AG', 'TC'] for pair in ref_query_base_pair]) > \
                min_edit - 1: # pass if there are AG or TC >= min_edit
            passCheck = True
            N_A2I = sum([pair in ['AG', 'TC'] for pair in ref_query_base_pair])
            return(passCheck, N_A2I)
    return(passCheck, 0)


def readBED(file):
    """
    file: should be \t delimited BED like file. 
    First 3 columns must be chr, 0-based start, end
    return a dataframe
    """
    import pandas as pd
    cols = ['chr', 'BEDstart', 'BEDend']
    df = pd.read_csv(file, sep="\t", names=cols, usecols=cols)
    return df


def main():
    import pysam as ps
    import argparse
    
    parser = argparse.ArgumentParser()
    # BED file from first pass A2G sites
    # e.g. /project2/yangili1/cdai/mRNA-editing/testPipe/Results/hs38/BCFCall/FinalAnno/ERR188060.bed
    parser.add_argument('--BED', type=str, required=True, 
                        help='BED file of the first pass editing sites') 
    # first STAR alignment, after recal
    # e.g. /project2/yangili1/cdai/mRNA-editing/testPipe/Results/hs38/BQSR/ERR188060_recal.bam
    parser.add_argument('--bamPre', type=str, required=True, 
                        help='bam file from first STAR alignment, likely after BQSR recal') 
    # second STAR alignment, after filtering
    # e.g. /project2/yangili1/cdai/mRNA-editing/testPipe/Results/hs38/Remap/ERR188060.filtered.bam 
    parser.add_argument('--bamPost', type=str, required=True,
                        help='bam file from second STAR alignment, after filtering') 
    parser.add_argument('--outTab', type=str, required=True,
                        help='out file name. Tab delimited')  
    parser.add_argument('--outBED', type=str, required=True,
                        help='output BED file, incl only chr, start, end')  
    args = parser.parse_args()



    # first 3 cols of BED, note pysam use 0-based thus directly using BED pos
    BED3 = readBED(args.BED)
    
    # get bam
    bam_pre = ps.AlignmentFile(args.bamPre, "rb")
    bam_pos = ps.AlignmentFile(args.bamPost, "rb")


    i = 0
    keep_idx, agreement, ad = [], [], [] # ad is A-G or T-C read depth
    for c, s in zip(BED3['chr'], BED3['BEDstart']):
        ck1, v1 = checkMapping(c, s, bam_pre, bam_pos)
        ck2, v2 = checkSite(c, s, bam_pos)

        if ck1 and ck2:
            keep_idx.append(i)
            agreement.append(v1)
            ad.append(v2)
        i += 1
    
    out_df = BED3.iloc[keep_idx]
    out_df.insert(3, 'agreement', agreement)
    out_df.insert(4, 'ad', ad)
    # write output
    out_df.to_csv(args.outTab, sep='\t', header=True, index=False)
    out_df.iloc[:, [0,1,2]].to_csv(args.outBED, sep='\t', header=False, index=False)

if __name__ == "__main__":
    main()
