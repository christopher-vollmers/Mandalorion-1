#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import argparse
import re
import os
import numpy as np
from os.path import isfile
import multiprocessing as mp
import mappy
import time



PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/utils/'
sys.path.append(os.path.abspath(PATH))

import SpliceDefineConsensus

parser = argparse.ArgumentParser()

parser.add_argument('--infile', '-i', type=str, action='store')
parser.add_argument('--path', '-p', type=str, action='store')
parser.add_argument('--cutoff', '-c', type=float, action='store')
parser.add_argument('--genome_file', '-g', type=str, action='store')
parser.add_argument('--splice_site_width', '-w', type=int, action='store')
parser.add_argument('--minimum_read_count', '-m', type=int, action='store')
parser.add_argument('--white_list_polyA', '-W', type=str, action='store')
parser.add_argument('--numThreads', '-n', type=str, action='store')
parser.add_argument('--junctions', '-j', type=str, action='store')
parser.add_argument('--upstream_buffer', '-u', type=str, action='store')
parser.add_argument('--downstream_buffer', '-d', type=str, action='store')
parser.add_argument('--abpoa', '-a', type=str, action='store')



if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()
infile = args.infile
out_path = args.path + '/'  # path where you want your output files to go
cutoff = float(args.cutoff)
genome_file = args.genome_file
splice_site_width = int(args.splice_site_width)
minimum_read_count = int(args.minimum_read_count)
white_list_polyA=args.white_list_polyA.split(',')
threads=int(args.numThreads)
junctions=args.junctions.split(',')
upstream_buffer=int(args.upstream_buffer)
downstream_buffer=int(args.downstream_buffer)
abpoa = args.abpoa


def process_locus(out_tmp,root,chrom,left_bounds_chrom, right_bounds_chrom,start,end,splice_site_width,minimum_read_count,junctions,cutoff,abpoa,verbose):
    infile=out_tmp+'/'+root+'.psl'
    peak_areas={}
    histo_left_bases, histo_right_bases, histo_cov,csDict = SpliceDefineConsensus.collect_reads(infile, chrom)
    if verbose:
        print(f'\t\tprocessing locus {chrom} {start} {end} covering {end-start:,} genomic bases and {len(csDict):,} sequencing reads' % (chrom,start,end)+' '*20)
#    else:
#        print(f'\t\tprocessing locus {chrom} {start} {end} covering {end-start:,} genomic bases and {len(csDict):,} sequencing reads' % (chrom,start,end)+' '*20, end='\r')
    peak_areas[chrom] = {}
    peak_areas[chrom]['l'] = {}
    peak_areas[chrom]['r'] = {}


    if verbose:
        print(f'\t\t\tcollecting annotated splice sites in locus {chrom} {start} {end}')
    peak_areas, toWrite_A_l = SpliceDefineConsensus.make_genome_bins(left_bounds_chrom, 'l', chrom, peak_areas,splice_site_width)
    peak_areas, toWrite_A_r = SpliceDefineConsensus.make_genome_bins(right_bounds_chrom, 'r', chrom, peak_areas,splice_site_width)
    if verbose:
        print(f'\t\t\tfinding high-confidence unannotated splice sites in locus {chrom} {start} {end}')
    peak_areas, toWrite_N_l = SpliceDefineConsensus.find_peaks(histo_left_bases[chrom], True, cutoff, histo_cov, 'l', peak_areas, chrom,csDict,start,end,splice_site_width,minimum_read_count,junctions)
    peak_areas, toWrite_N_r = SpliceDefineConsensus.find_peaks(histo_right_bases[chrom], False, cutoff, histo_cov, 'r', peak_areas, chrom,csDict,start,end,splice_site_width,minimum_read_count,junctions)

    peakCounter={}
    peakCounter['l']=0
    peakCounter['r']=0
    spliceDict={}
    spliceDict[chrom]={}
    for toWrite in [toWrite_A_l,toWrite_A_r,toWrite_N_l,toWrite_N_r]:
       for chrom,start,end,type1,side,prop in toWrite:
            peakCounter[side]+=1
            peaks=str(peakCounter[side])
            splice_left = int(start)
            splice_right = int(end)
            for base in np.arange(splice_left, splice_right + 1):
                spliceDict[chrom][base] = type1 + side + peaks

    if verbose:
        print(f'\t\t\tsorting read alignments into splice junctions for locus {chrom} {start} {end}')
    start_end_dict, start_end_dict_mono = SpliceDefineConsensus.sort_reads_into_splice_junctions(spliceDict, infile)

    if verbose:
        print(f'\t\t\tfinding TSS and polyA sites for splice junctions chains for locus {chrom} {start} {end}')
    seqDict=SpliceDefineConsensus.define_start_end_sites(start_end_dict, start_end_dict_mono,upstream_buffer,downstream_buffer,minimum_read_count)
    IsoData={}
    if verbose:
        print(f'\t\t\tcreating consensus sequences for {len(seqDict):,} putative, unfiltered isoforms for locus {chrom} {start} {end}')
    for isoform,reads in seqDict.items():
        consensus,names=SpliceDefineConsensus.determine_consensus(reads,out_tmp+'/'+root,abpoa)
        IsoData[isoform]=[consensus,names]
    if verbose:
        print(f'\t\t\tfinished processing of locus {chrom} {start} {end}')


    return IsoData

def main():

    out_tmp = out_path + '/tmp_SS'

    left_bounds, right_bounds = {}, {}
    if genome_file!='None':
        if genome_file.endswith('.gtf.gz') or genome_file.endswith('.gtf'):
            print('\tparsing annotated splice sites')
            chrom_list, left_bounds, right_bounds, polyAWhiteList = SpliceDefineConsensus.parse_genome(genome_file, left_bounds, right_bounds, white_list_polyA)
        else:
            print('\tgenome annotation file does not end with .gtf or .gtf.gz. File will be ignored. Splice sites will be entirely read derived and no polyA sites will be white-listed')
            chrom_list=set()
            polyAWhiteList=[]
    else:
        print('\tNo genome annotation provided, so splice sites will be entirely read derived and no polyA sites will be white-listed')
        chrom_list=set()
        polyAWhiteList=[]

    outPolyA=open(out_path+'/polyAWhiteList.bed','w')

    if '0' not in white_list_polyA:
        print('\t'+str(len(polyAWhiteList)), 'poly(A) sites whitelisted')
        for chrom,direction,end,transcript_id in polyAWhiteList:
            polyA=int(end)
            polyAstart=polyA-20
            polyAend=polyA+20
            outPolyA.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,str(polyAstart),str(polyAend),transcript_id,'0',direction))
    outPolyA.close()

    print('\tcollecting loci')
    chrom_list,roots = SpliceDefineConsensus.get_parsed_files(out_tmp,chrom_list)
    numThreads=threads
    roots=sorted(list(roots),key=lambda x: (x.split('~')[0],int(x.split('~')[1])))
    out=open(out_path+'/Isoform_Consensi.fasta','w')
    out_r2i=open(out_path+'/reads2isoforms.txt','w')
    results={}
    print('\tstarting multithreaded processing of loci')
    pool = mp.Pool(numThreads,maxtasksperchild=1)
    for root in roots:
        chrom,start,end=root.split('~')
        start=int(start)
        end=int(end)
        left_bounds_sub,right_bounds_sub=SpliceDefineConsensus.prepare_locus(chrom,start,end,left_bounds,right_bounds)
        results[root]=pool.apply_async(process_locus,[out_tmp,root,chrom,left_bounds_sub[chrom], right_bounds_sub[chrom],start,end,splice_site_width,minimum_read_count,junctions,cutoff,abpoa,False])

    pool.close()
    previous=0
    delay=True
    while delay:
        time.sleep(60)
        total_roots=0
        finished_roots=0
        for root in results:
            total_roots+=1
            if results[root].ready():
                finished_roots+=1
        if previous==finished_roots:
            delay=False
        previous=finished_roots
        print(f'\t\tfinished {finished_roots} of {total_roots} loci',' '*30)


    print(f'\n\t{total_roots-finisihed_roots} loci took too long to complete\n\tterminating multi-threading pool',' '*30)
    pool.terminate()
    pool.join()

    counter=0
    print('\twriting isoform sequences to file for loci that completed')
    unfinished_roots=set()
    for root in roots:
        chrom,start,end=root.split('~')
        if results[root].ready():
            IsoData = results[root].get()
            for isoform in IsoData:
                counter+=1
                consensus,names = IsoData[isoform]
                nameString='Isoform'+str(counter)+'_'+str(len(names))
                out.write('>%s\n%s\n' % (nameString,consensus))
                for name in names:
                    out_r2i.write('%s\t%s\n' %(name,nameString))
        else:
            unfinished_roots.add(root)

    print('\tperforming single-threaded processing for loci that did not complete')
    results={}
    unfinished_roots=sorted(list(unfinished_roots),key=lambda x: (x.split('~')[0],int(x.split('~')[1])))
    for root in unfinished_roots:
        print('processing',root, 'single-threaded')
        chrom,start,end=root.split('~')
        start=int(start)
        end=int(end)
        left_bounds_sub,right_bounds_sub=SpliceDefineConsensus.prepare_locus(chrom,start,end,left_bounds,right_bounds)
        results[root]=process_locus(out_tmp,root,chrom,left_bounds_sub[chrom], right_bounds_sub[chrom],start,end,splice_site_width,minimum_read_count,junctions,cutoff,abpoa,True)
    print('\twriting isoform sequences to file for single-threaded loci')
    for root in unfinished_roots:
        chrom,start,end=root.split('~')
        IsoData = results[root]
        for isoform in IsoData:
            counter+=1
            consensus,names = IsoData[isoform]
            nameString='Isoform'+str(counter)+'_'+str(len(names))
            out.write('>%s\n%s\n' % (nameString,consensus))
            for name in names:
                out_r2i.write('%s\t%s\n' %(name,nameString))




    out.close()
    out_r2i.close()

main()
