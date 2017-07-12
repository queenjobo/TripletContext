'''
@queenjobo
July 2017

This script takes a tab file with first 6 columns ID BAMFILE CHR POS REF ALT and gets
the k (default triplet) context around the variant and output that context and the proportion
of reads that agree with that context as a quality metric

INPUT FILE: tab file with first 6 columns ID BAMFILE CHR POS REF ALT 
OUTPUT FILE: tab file with same as input file with added columns CONTEXT CONTEXT_QUALITY

NOTE: can accept bam or cram files

USAGE: python get_triplet_context.py --variants VARIANT_FILE.tab [--k K] [--output OUTPUT_FILE]
EXAMPLE: python get_triplet_context.py --variants variants_example.tab

python get_triplet_context.py --variants /lustre/scratch115/projects/ddd/users/jk18/mutatpheno/DNG/DNG.output.filtered.29.06.17_forIGV.tab --output /lustre/scratch115/projects/ddd/users/jk18/mutatpheno/DNG/DNG.output.filtered.29.06.17_context.tab

'''

#-----------------IMPORTS-------------------

import os
import argparse
import pysam
from collections import Counter

#---------------FUNCTIONS-------------------

def get_options():
    '''parse command line options'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants',type = str,help = 'path to tab file of variants with first 6 columns ID BAMFILE CHR POS REF ALT',required = True)
    parser.add_argument('--k', type = int, help = "size of context to retrieve (bp)", required = False, default = 3)
    parser.add_argument('--output', type = str, default = "myoutput.tab", help = "output file path", required = False)
    args = parser.parse_args()
    return(args)

def is_cram(myfile):
    out = None
    if myfile.endswith(".cram"):
        out = True
    elif myfile.endswith(".bam"):
        out = False
    return(out)
    
def get_kmer_read(read,pos,kside):
    '''given read from pysam, pos and kside (number around kmer) get kmer string'''
    myseq = ""
    start = (pos-read.reference_start)-kside-1
    stop = (pos-read.reference_start)+kside
    if start > 0 and stop > 0:
        myseq = read.query_sequence[start:stop]
    return(myseq)

def get_most_common(mylist):
    mycount = Counter(mylist)
    return(mycount.most_common(1))

def get_context(k,bamfile,chrm,pos,ref,alt):
    '''for variant given bamfile and variant details get context'''
    kmers = []
    kside = int((k-1)/2)
    #open cram/bam file
    if is_cram(bamfile):
        samfile = pysam.AlignmentFile(bamfile,'rc')
    else:
        samfile = pysam.AlignmentFile(bamfile,'rb')
    iter = samfile.fetch(str(chrm),pos-kside,pos+kside)
    for read in iter:
        seq = get_kmer_read(read,pos,kside)
        #only if has alt
        if len(seq) == k: 
            if seq[kside] == alt:
                kmers.append(seq)
    totkmer = len(kmers) 
    if totkmer == 0:
        altcontext = "NA"
        refcontext ="NA"
        seqcontext = "NA"
        prop ="NA"
    else:
        mostcommon = get_most_common(kmers)
        altcontext = mostcommon[0][0]
        prop = float(mostcommon[0][1])/len(kmers)
        refcontext = get_ref_context(altcontext,kside,ref)
    return(refcontext,altcontext,prop)

def get_ref_context(altcontext,kside,ref):
    refcontext = list(altcontext)
    refcontext[kside] = ref
    return("".join(refcontext))

def get_all_contexts(variants,k):
    '''go through all variants and get context'''
    newlines = []
    with open(variants,'r') as f:
        for line in f:
            #extract details from line
            fields = line.strip("\n").split("\t")
            bamfile = fields[1].split(",")[0]
            chrm = fields[2]
            pos = int(fields[3])
            ref = fields[4]
            alt = fields[5]
            #get context for variant
            refcontext,altcontext,prop = get_context(k,bamfile,chrm,pos,ref,alt)
            newline = "\t".join(fields + [refcontext,altcontext,str(prop)])
            newlines.append(newline)
    return(newlines)
            
def write_to_file(mylines,output):
    with open(output,'w') as nf:
        nf.write("\n".join(mylines)+"\n")

def main():
    #parse options
    args = get_options()
    alllines = get_all_contexts(args.variants,args.k)
    write_to_file(alllines,args.output)
    
#---------------SCRIPT----------------------
if __name__ == "__main__":
    main()