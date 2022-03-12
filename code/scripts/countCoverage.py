'''
@author: Chao Dai
data: 1/24/2022
Description: provide a bed file of single positions, compute coverages
'''



def countCoverage(bam, ch, st):
    '''
    bam: pysam AlignmentFile object
    ch: chr
    st: start position, 0 based
    en: end position 0 based
    Description: given a site and a sample bam, return a dictionary describing coverage of each base A, C, G, T
    '''
    coverage = bam.count_coverage(ch, st, st+1, quality_threshold=20)
    A, C, G, T = coverage
    return {'A': A[0], 'C':C[0], 'G':G[0], 'T':T[0]}


def combineCoverage(bamFile, bedFile):
    '''
    bam: bamFile path
    bed: 3 column bed file format
    return a dataframe with the first 3 columns identical to the bedFile
    but append the count A, C, G, T at each position

    return: a dataframe with same rows as the original bed, but 
    adding 4 columns of A, C, G, T, coverage
    '''
    import pandas as pd
    import pysam as ps
    from datetime import datetime
    print(f'## {datetime.now()}, counting coverage...')

    # read in bed
    bed = pd.read_csv(bedFile, sep="\t", 
                      names=['chr', 'BEDstart', 'BEDend'], 
                      usecols=[0,1,2], 
                      dtype={'chr': str, 'BEDstart': int, 'BEDend': int})
    
    # construct pysam object then compute coverage
    with ps.AlignmentFile(bamFile, 'r') as bam:
        cov = list(map(lambda ch, st: countCoverage(bam, ch, st), bed.chr, bed.BEDstart))
    
    cov = pd.DataFrame(cov) # covert dict to df
    
    print(f'## {datetime.now()} done.')
    return pd.concat([bed, cov], axis=1)

def main():

    import argparse
    import pandas as pd

    # set up argument parser object
    parser = argparse.ArgumentParser()
    parser.add_argument('--BAM', type=str, required=True,
                        help='bam file from which to count reads') 
    parser.add_argument('--BED', type=str, required=True,
                        help='BED of filtered A>I sites from the last step, e.g. removeHomopolymers') 
    parser.add_argument('--outTab', type=str, required=True,
                        help='out file name. Tab delimited')  
    args = parser.parse_args()

    coverage = combineCoverage(args.BAM, args.BED)
    coverage.to_csv(args.outTab, sep='\t', index=False)


if __name__ == "__main__":
    main()