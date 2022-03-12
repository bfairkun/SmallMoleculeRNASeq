'''
author: Chao Dai
date: 1/10/22
Description: check 5' reads and remove sites that are in 5' reads
'''


def readBED(file):
    """
    file: should be \t delimited BED like file. First 3 columns must be chr, 0-based start, end
    return a dataframe
    """
    import pandas as pd
    cols = ['chr', 'BEDstart', 'BEDend']
    df = pd.read_csv(file, sep="\t", names=cols, usecols=cols)
    return df


def check5prime(ch, s, bam):
    '''
    ch: str datatype, e.g. 'chr1'
    s: int datatype, e.g. 68455, BED format, 0-based inclusive
    bam: ps.AlignmentFile
    description: get all reads that support the mismatch in 6bp of 5'
    returns a dictionary of lists.
    '''
    e = s + 1 # end pos, 0-based open
    
    # for each site, find the reads that support a mismatch AND that the mismatch
    # is within 6bp of 5'
    failed_reads = {'query_name':[], 'is_read1':[], 'mismatch_qpos':[]}
    for read in bam.fetch(ch, s, e):
        aligned_pair = [tup for tup in read.get_aligned_pairs(matches_only=True, with_seq=True) if tup[1] == s]
    try:
        query_base_pos, ref_base_pos, ref_base = aligned_pair[0]
        if ref_base.upper() != ref_base:
            if (read.is_read1 and query_base_pos < 6) or (read.is_read2 and read.query_length - query_base_pos < 6): # also check last 6bp if read2
                failed_reads.get('query_name').append(read.query_name)
                failed_reads.get('is_read1').append(read.is_read1)
                failed_reads.get('mismatch_qpos').append(query_base_pos)
    except IndexError:
        pass
    
    return failed_reads


def Calc5PrimePercent(bed, bam): 
    '''
    bed: pd.Dataframe of Bed format file, 3 columns
    bam: ps.AlignmentFile object
    description: compute the percentage of reads that overlap a site within 6bp of 5'
    return value: a list of percentages. Each number cooresponds a given site (row in the bed dataframe)
    '''
    import checkRemap as remap
    
    c = [] # store percentage
    for ch, st in zip(bed.chr, bed.BEDstart):
        a = len(check5prime(ch, st, bam).get('query_name')) # get all reads that support the mismatch in 6bp of 5', compute how many of such reads
        b = remap.checkSite(ch, st, bam)[1] # number of reads supporting the mismatch
        c.append(a/b) # ratio: out of all reads supporting the mismatch, how many are in the 6bp of 5'
    return c


def main():
    import pysam as ps
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--BED', type=str, required=True, help='BED file of editing sites, requires the first 3 columns') 
    # e.g. /project2/yangili1/cdai/mRNA-editing/testPipe/Results/hs38/Remap/ERR188060.filtered.bam 
    parser.add_argument('--bam', type=str, required=True, help='bam file from second STAR alignment, after filtering') 
    parser.add_argument('--th', type=float, required=False, default=0.2, help="maximum threshold of percent of reads supporting site within first bp of 5'. Default 0.2")
    parser.add_argument('--outBED', type=str, required=True, help='output BED file, incl only chr, start, end')  
    args = parser.parse_args()

    BED = readBED(args.BED)
    BAM = ps.AlignmentFile(args.bam, "rb")
    THRESHOLD = args.th

    # compute ratio of reads wtihin 6bp of 5' for each site
    PERC = Calc5PrimePercent(BED, BAM)
    PASS_INDX = [p < THRESHOLD for p in PERC]

    # write passing sites to a BED file
    BED[PASS_INDX].to_csv(args.outBED, sep='\t', header=False, index=False)



if __name__ == '__main__':
    main()