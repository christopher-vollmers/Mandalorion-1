#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import argparse
import re
import os
import numpy as np
import mappy as mp
import gzip
import time

def prepare_locus(chrom,start,end,left_bounds,right_bounds):
    if chrom not in left_bounds:
        left_bounds[chrom] = {'5': [], '3': []}
    if chrom not in right_bounds:
        right_bounds[chrom] = {'5': [], '3': []}

    left_bounds_sub = {}
    right_bounds_sub = {}
    left_bounds_sub[chrom] =  {'5': [], '3': []}
    right_bounds_sub[chrom] = {'5': [], '3': []}
    for side in ['5','3']:
        for pos in left_bounds[chrom][side]:
            if start<pos<end:
                left_bounds_sub[chrom][side].append(pos)
        for pos in right_bounds[chrom][side]:
            if start<pos<end:
                 right_bounds_sub[chrom][side].append(pos)
    return left_bounds_sub,right_bounds_sub



def clean_psl(psl_file, clean_psl_file,primary):

    minimum_intron_size=10

    used=set()
    out = open(clean_psl_file, 'w')
    for line in open(psl_file):
        a = line.strip().split('\t')
        start = int(a[15])
        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]

        name=a[9]
        if primary and name in used:
            continue

        blockstarts_clean = []
        readstarts_clean = []
        blocksizes_clean = []
        size_gap = []

        for x in range(len(blocksizes)):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            blockend = blockstart + blocksize
            size_gap.append(blocksize)

            try:
                next_blockstart = int(blockstarts[x + 1])
                gap = next_blockstart - blockend
                size_gap.append(gap)
            except IndexError:
                pass

        new_size_gap = []
        block = 0
        for index in range(len(size_gap)):
            if index % 2 == 0:
                block += size_gap[index]
            if index % 2 == 1:
                if size_gap[index] < minimum_intron_size:
                    block += size_gap[index]
                else:
                    new_size_gap.append(block)
                    new_size_gap.append(size_gap[index])
                    block = 0

        new_size_gap.append(block)

        current_position = start
        current_readposition = int(a[11])
        for index in range(0, len(new_size_gap), 1):
            inter = new_size_gap[index]
            if index % 2 == 0:
                blockstarts_clean.append(str(current_position))
                blocksizes_clean.append(str(inter))
                readstarts_clean.append(str(current_readposition))

                current_position += inter
                current_readposition += inter

            if index % 2 == 1:
                current_position += inter

        a[17] = str(len(blockstarts_clean))

        blockstarts_clean.append('')
        blocksizes_clean.append('')
        readstarts_clean.append('')

        a[18] = (',').join(blocksizes_clean)
        a[19] = (',').join(readstarts_clean)
        a[20] = (',').join(blockstarts_clean)



        new_line = ('\t').join(a)
        out.write(new_line + '\n')
        used.add(name)



def get_parsed_files(out_tmp,chrom_list):
    roots=set()
    for file1 in os.listdir(out_tmp):
        if os.path.isfile(out_tmp+'/'+file1):
            if '.psl' in file1:
                root=file1.split('.psl')[0]
                chrom,start,end = root.split('~')
                chrom_list.add(chrom)
                roots.add(root)
    return chrom_list,roots

def getCSaroundSS(CS,genomePosition,start,end):
    '''
    Goes through the CIGAR string by pairing numbers with the letters.
    Counts the number of matches/mismatches and indels.
    Returns the total number of all, matches/mismatches (M), and indels.
    '''

    status=''
    p = re.compile(r'([\\=\\+\-\\*\\~])')
    splitCstr = [i+j for i,j in zip(p.split(CS)[1::2], p.split(CS)[2::2])]
    record=[]
    spliceIndex=False
    for i in range(len(splitCstr)):
        sub=splitCstr[i]
        status=sub[0]
        entry=sub[1:]
        if status == '=':
            for item in entry:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '+':
            for item in entry:
                record.append((status,item))
        if status == '-':
            for item in entry:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '*':
            for item in entry[::2]:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '~':
            intronLength=int(entry[2:-2])
            genomePosition+=intronLength
            record.append('|'+entry+'|')
            if start<=genomePosition<=end:
                  spliceIndex=len(record)
    bases='nnnn'
    left=''
    right=''

    if spliceIndex:
        for index in range(max(spliceIndex-10,0),min(spliceIndex+10,len(record))):
            item=record[index]
            if '|' in item:
               bases=item[1:3]+item[-3:-1]
               left = record[index-5:index]
               right = record[index+1:index+6]
    return bases,left,right

def scan_for_best_bin(entry, dist_range, density_dict, peak_areas, chrom, side):
    best_extra_list, peak_center = 0, 0
    cov_area = {}
    best_dirs = {'+': 0, '-': 0}
    best_names = []
    for x in dist_range:
        extra_list_exp = 0
        dirs = {'+': 0, '-': 0}
        names=[]
        cov_set = {}
        called = False
        for y in dist_range:
            if entry + x + y in peak_areas[chrom][side]:
                called = True
        if called:
            continue
        for y in dist_range:
            if entry + x + y in density_dict:
                for item in density_dict[entry + x + y]:
                    extra_list_exp += 1
                    dirs[item[4]] += 1
                    names.append(item[0])
                    for covered_position in item[3]:
                        if covered_position not in cov_set:
                            cov_set[covered_position] = 0
                        cov_set[covered_position] += 1

        if extra_list_exp > best_extra_list:
            best_extra_list = extra_list_exp
            best_names=names
            peak_center = entry + x
            cov_area = cov_set
            best_dirs = dirs

    return best_extra_list, peak_center, cov_area, best_dirs,best_names


def determine_cov(cov_area, chrom, reverse, peak_center, histo_cov):
    cov = [0]
    cov_area2 = []
    for covered_position, count in cov_area.items():
        if count > 1:
            cov_area2.append(covered_position)

    cov_area = sorted(cov_area2, reverse=reverse)
    counter = 0
    for base_f in cov_area:
        count = False
        if reverse and base_f < peak_center:
            count = True
        elif not reverse and base_f > peak_center:
            count = True
        if count:
            if counter <= 3:
                counter += 1
                base_f = myround(base_f)
                if base_f in histo_cov[chrom]:
                    cov.append(histo_cov[chrom][base_f])
            else:
                break
    cov = max(cov)
    return cov, cov_area


def myround(x, base=10):
    '''Rounds to the nearest base'''
    return int(base * round(float(x) / base))


def find_peaks(density_dict, reverse, cutoff, histo_cov, side, peak_areas, chrom, csDict,start,end,splice_site_width,minimum_read_count,junctions):
    dist_range=[0]
    toWrite=[]
    for shift in range(1,splice_site_width+1):
        dist_range.append(shift)
        dist_range.append(-shift)

    entry_list = []
    for entry,density in density_dict.items():
        if len(density)>=minimum_read_count:
            entry_list.append([entry, density])
    for entry, density in sorted(entry_list, key=lambda x: len(x[1]), reverse=True):
        if entry in peak_areas[chrom][side]:
            continue
        best_extra_list, peak_center, cov_area, best_dirs,best_names = scan_for_best_bin(entry, dist_range, density_dict, peak_areas, chrom, side)
        cov, cov_area = determine_cov(cov_area, chrom, reverse, peak_center, histo_cov)

        if cov > 0:
            proportion = round(best_extra_list / cov, 3)
            if proportion > cutoff:
                Left_to_Right = best_dirs['+']
                Right_to_Left = best_dirs['-']
                Type = ''

                if Left_to_Right < Right_to_Left:
                    Type = '3' if reverse else '5'

                elif Left_to_Right > Right_to_Left:
                    Type = '5' if reverse else '3'

                if not Type:
                    continue
                passed = characterize_splicing_event([chrom, peak_center - splice_site_width, peak_center + splice_site_width], best_names, csDict,junctions)
                if passed:

                    start=str(peak_center - splice_site_width)
                    end=str(peak_center + splice_site_width)
                    toWrite.append([chrom,start,end,Type,side,str(proportion)])

                    for base in range(int(start), int(end) + 1):
                        peak_areas[chrom][side][base] = 1
        else:
            continue
    return peak_areas, toWrite


def collect_reads(reads, target_chrom):
    histo_cov, histo_left_bases, histo_right_bases, subsample = {}, {}, {}, 100000
    histo_cov[target_chrom] = {}
    histo_left_bases[target_chrom] = {}
    histo_right_bases[target_chrom] = {}
    csDict={}
    lineCounter=0
    for line in open(reads):
        lineCounter+=1


    indexes = np.random.choice(range(0, lineCounter, 1),
                               size=min(lineCounter,subsample),
                               replace=False)

    indexes=set(indexes)

    lineCounter=-1
    for line in open(reads):
        lineCounter+=1
        if lineCounter in indexes:
            a = line.strip().split('\t')
            chrom = a[13]
            if chrom==target_chrom:
                name = a[9]
                dirn = a[8]
                length = int(a[10])
                begin, span = int(a[15]), int(a[16])
                blocksizes = a[18].split(',')[:-1]
                blockstarts = a[20].split(',')[:-1]
                accuracy = float(a[21])
                cs=a[22]
                csDict[name]=(cs,begin)
                cov_set = set()
                low_bounds, up_bounds = [], []
                aligned_bases = 0
                for x in range(0, len(blocksizes)):
                    blockstart = int(blockstarts[x])
                    blocksize = int(blocksizes[x])
                    aligned_bases += blocksize
                    blockend = blockstart + blocksize
                    for y in range(0, blocksize, 10):
                        rounded = myround(blockstart + y)
                        cov_set.add(rounded)
                    for yy in range(y, blocksize):
                        rounded = myround(blockstart + yy)
                        cov_set.add(rounded)
                    if blockstart != begin:
                        up_bounds.append(blockstart)
                    if blockend != span:
                        low_bounds.append(blockend)

                for rounded in cov_set:
                    if rounded not in histo_cov[chrom]:
                        histo_cov[chrom][rounded] = 0
                    histo_cov[chrom][rounded] += 1

                if accuracy<0.9:
                    continue
                for low_bound in low_bounds:
                    if low_bound not in histo_left_bases[chrom]:
                        histo_left_bases[chrom][low_bound] = []
                    histo_left_bases[chrom][low_bound].append([name, begin, span, cov_set, dirn,accuracy,round(aligned_bases,-2)])
                for up_bound in up_bounds:
                    if up_bound not in histo_right_bases[chrom]:
                        histo_right_bases[chrom][up_bound] = []
                    histo_right_bases[chrom][up_bound].append([name, begin, span, cov_set, dirn,accuracy,round(aligned_bases,-2)])
    return histo_left_bases, histo_right_bases, histo_cov,csDict


def parse_genome(input_file, left_bounds, right_bounds, white_list_polyA):
    polyAWhiteList=[]
    chrom_list = set()
    gene_dict = {}
    if input_file.endswith('.gtf.gz'):
        print('\t\tgtf file ends on .gz and will be treated as gzipped')
        input=gzip.open(input_file,'rt')
    elif input_file.endswith('.gtf'):
        input=open(input_file,'r')

    for line in input:
        pAwl=False
        for element in white_list_polyA:
            if element in line:
                pAwl=True
        a = line.strip().split('\t')
        if len(a) <= 7:
            continue
        if a[2] == 'exon':
            testKey = a[8].split('transcript_id "')[1].split('"')[0]
            if not gene_dict.get(testKey):
                gene_dict[testKey] = []
            gene_dict[testKey].append((a[0], a[3], a[4], a[6], pAwl))

    for transcript_id in gene_dict:
        transcript_data = gene_dict[transcript_id]

        chrom = transcript_data[0][0]
        direction=transcript_data[0][3]
        chrom_list.add(chrom)
        if chrom not in right_bounds:
            left_bounds[chrom] = {'5': [], '3': []}
            right_bounds[chrom] = {'5': [], '3': []}

        start = sorted(transcript_data, key=lambda x: int(x[1]))[0][1]
        end = sorted(transcript_data, key=lambda x: int(x[2]), reverse=True)[0][2]
        direction=transcript_data[0][3]
        pAwl=transcript_data[0][4]
        if pAwl:
            if direction=='+':
                polyAWhiteList.append((chrom,direction,end,transcript_id))
            elif direction=='-':
                polyAWhiteList.append((chrom,direction,start,transcript_id))

        for entry in transcript_data:
            if entry[1] != start:
                if entry[3] == '+':
                    right_bounds[chrom]['3'].append(int(entry[1]) - 1)
                elif entry[3] == '-':
                    right_bounds[chrom]['5'].append(int(entry[1]) - 1)
            if entry[2] != end:
                if entry[3] == '+':
                    left_bounds[chrom]['5'].append(int(entry[2]))
                if entry[3] == '-':
                    left_bounds[chrom]['3'].append(int(entry[2]))
    return chrom_list, left_bounds, right_bounds, polyAWhiteList


def make_genome_bins(bounds, side, chrom, peak_areas,splice_site_width):
    toWrite=[]
    for type1 in ['5', '3']:
        covered = {}
        position_list = sorted(bounds[type1], key=int)
        for index1 in range(0, len(position_list)):
            if index1 in covered:
                continue
            sub_list = [position_list[index1]]
            for index2 in range(index1, len(position_list)):
                if position_list[index2] - max(sub_list) <= splice_site_width:
                    sub_list.append(position_list[index2])
                    covered[index2] = 1
                else:
                    break
            single = False
            if len(sub_list) > 1:
                splice_dists = []
                for splice_pos in range(0, len(sub_list) - 1):
                    splice_dists.append(int(sub_list[splice_pos + 1]) - int(sub_list[splice_pos]))
                if min(splice_dists) > 3:
                    for x in range(0, len(sub_list), 1):
                        if x != 0:
                            start = int(sub_list[x] - ((sub_list[x] - sub_list[x - 1]) / 2))
                        else:
                            start = int(sub_list[x]) - splice_site_width
                        if x != len(sub_list) - 1:
                            end = int(sub_list[x] + ((sub_list[x + 1] - sub_list[x]) / 2))
                        else:
                            end = int(sub_list[x]) + splice_site_width

                        toWrite.append([chrom,str(start),str(end),type1,side,'A'])
                        for base in range(start, end + 1):
                            peak_areas[chrom][side][base] = 1

                else:
                    single = True
            else:
                single = True
            if single:
                start = min(sub_list) - splice_site_width
                end = max(sub_list) + splice_site_width
                toWrite.append([chrom,str(start),str(end),type1,side,'A'])
                for base in range(start, end + 1):
                    peak_areas[chrom][side][base] = 1

    return peak_areas,toWrite


def prune_locus(previous_start,previous_end,reads,minCount):
    minCount=int(minCount)
    coords=[]
    for line in reads:
        a=line.strip().split('\t')
        chrom=a[13]
        start=int(a[15])
        end=int(a[16])
        name=a[9]
        coords.append((start,1,line))
        coords.append((end,-1,line))
    coords=sorted(coords)
    coverage=0

    lines=[]
    locus=False
    loci=[]
    for coord in coords:
         loc,value,line=coord[0],coord[1],coord[2]
         coverage+=value
         if locus:
             if coverage<minCount:
                 locusEnd=loc
                 root=f'{chrom}~{locusStart}~{locusEnd}'
                 rootReads=lines
                 loci.append((root,rootReads))
                 lines=[]
                 locus=False
         else:
             if coverage >= minCount:
                 locus=True
                 locusStart=loc
         if value == 1:
             lines.append(line)
    return loci


def get_loci(infile,out_tmp,minCount):
    chrom_list=set()
    outDict={}
    reads=[]
    new=False
    previous_chrom=''
    previous_start=0
    previous_end=0
    roots=set()
    total_psl=0
    print('\t\tparsing psl and splitting into loci')
    for line in open(infile):
        a=line.strip().split('\t')
        chrom=a[13]
        chrom_list.add(chrom)
        start=int(a[15])
        end=int(a[16])
        if chrom!=previous_chrom or start>previous_end:
            new=True

        if not new:
            previous_end=max(end,previous_end)
            reads.append(line)
        else:
            if reads:
                root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
                loci=prune_locus(previous_start,previous_end,reads,minCount)
                print(f'\t\tPreliminary locus {root} with {len(reads)} reads split into {len(loci)} final loci',' '*20,end='\r')
                for root,root_reads in loci:
#                    print(f'\t\t\t{root} {len(root_reads)}')
                    fh=open(out_tmp+'/'+root+'.psl', 'w')
                    roots.add(root)

                    for read in root_reads:
                        fh.write(read)
                        total_psl+=1
                        name=read.strip().split('\t')[9]
                        outDict[name]=root
                    fh.close()
                reads=[]
            reads.append(line)   ### bug in v4.0.0. Left the first read in each locus unused
            previous_chrom=chrom
            previous_end=end
            previous_start=start
            new=False
    if reads:
        loci=prune_locus(previous_start,previous_end,reads,minCount)
        for root,root_reads in loci:
            fh=open(out_tmp+'/'+root+'.psl', 'w')
            roots.add(root)
            for line in reads:
                total_psl+=1
                fh.write(line)
                name=line.strip().split('\t')[9]
                outDict[name]=root
            fh.close()

    print(f'\t\tsplit {total_psl} psl entries into {len(roots)} loci',' '*20)



def get_chromosomes(infile,out_tmp,fastaList):
    chrom_list=set()
    outDict={}
    reads=[]
    new=False
    previous_chrom=''
    previous_start=0
    previous_end=0
    roots=set()
    total_psl=0
    print('\t\tparsing psl and splitting into loci')
    for line in open(infile):
        a=line.strip().split('\t')
        chrom=a[13]
        chrom_list.add(chrom)
        start=int(a[15])
        end=int(a[16])
        if chrom!=previous_chrom or start>previous_end:
            new=True

        if not new:
            previous_end=max(end,previous_end)
            reads.append(line)
        else:
            if reads:
                root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
                print('\t\tdefining locus',(' ').join(root.split('~')),' '*20,end='\r')

                fh=open(out_tmp+'/'+root+'.psl', 'w')
                roots.add(root)

                for read in reads:
                    fh.write(read)
                    total_psl+=1
                    name=read.strip().split('\t')[9]
                    outDict[name]=root
                fh.close()
                reads=[]
            reads.append(line)   ### bug in v4.0.0. Left the first read in each locus unused
            previous_chrom=chrom
            previous_end=end
            previous_start=start
            new=False
    if reads:
        root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
        fh=open(out_tmp+'/'+root+'.psl', 'w')
        roots.add(root)
        for line in reads:
            total_psl+=1
            fh.write(line)
            name=line.strip().split('\t')[9]
            outDict[name]=root
        fh.close()

    print('\t\tsplit', total_psl, 'psl entries into loci',' '*20)



def characterize_splicing_event(a,names,csDict,junctions):
    passed = False
    chromosome = a[0]
    splice_left = int(a[1])
    splice_right = int(a[2])
    CSs=[]
    indeces = np.random.choice(np.arange(0, len(names)),
                                   min(len(names), 500), replace=False)
    for index in indeces:
        name = names[index]
        CSs.append(csDict[name])
    basecontext={}
    basecontext['allowed']=0
    basecontext['all']=0
    leftCS={}
    leftCS['Total']=0
    rightCS={}
    rightCS['Total']=0
    for status in ['*','+','-','=','|']:
        leftCS[status]=0
        rightCS[status]=0
    for cs,genomePosition in CSs:
        bases,left,right = getCSaroundSS(cs,genomePosition,splice_left,splice_right)
        basecontext['all'] += 1
#        if bases == 'gtag' or bases == 'ctac':
        if bases in junctions:
            basecontext['allowed'] += 1
        else:
            if bases not in basecontext:
                basecontext[bases]=0
            basecontext[bases]+=1
        if left:
            for thing in left:
                leftstatus = thing[0]
                leftCS[leftstatus]+=1
                leftCS['Total']+=1
        if right:
            for thing in right:
                rightstatus = thing[0]
                rightCS[rightstatus]+=1
                rightCS['Total']+=1

    basepass = False
    if basecontext['allowed']/basecontext['all']>0.85:
        basepass = True
        leftacc=leftCS['=']/leftCS['Total']
        rightacc=rightCS['=']/rightCS['Total']
    if basepass:
        if leftacc>0.85 and rightacc>0.85:
            passed = True
    CIGARs, CigarDict = 0, 0
    return passed



def find_ends(starts, ends, identity, count_dict, upstream_buffer,downstream_buffer,minimum_feature_count):
    start_peaks, end_peaks, extend = {}, {}, True
    if '+' in identity:
        count_dict['+'].add(identity)
    else:
        count_dict['-'].add(identity)


    start_count = {}
    end_count = {}
    for position in starts:
        if position not in start_count:
            start_count[position] = 0
        start_count[position] += 1

    for position in ends:
        if position not in end_count:
            end_count[position] = 0
        end_count[position] += 1

    show=False
    for position in sorted(starts):
        if position - upstream_buffer not in start_peaks:
            window_count = 0
            for i in np.arange(0,10):
                window_position = position + i
                if window_position in start_count:
                    window_count += start_count[window_position]
            if window_count >= minimum_feature_count:
                original_bin=[]
                for shift in np.arange(-upstream_buffer, downstream_buffer):
                    start_peaks[position + shift] = position
                    original_bin.append(position + shift)
                if extend:
                    bins=[]
                    for i in np.arange(min(original_bin),max(original_bin)):
                        bin=0
                        for shift in np.arange(0,10):
                            if i+shift in start_count:
                                 bin+=start_count[i+shift]
                        bins.append(bin)

                    best_bin=max(bins)
                    extended=True
                    adjacent=position - upstream_buffer
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent - i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in start_count:
                                window_count += start_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in start_peaks:
                                    start_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['start_left'].add(identity)

                    extended = True
                    adjacent = position + downstream_buffer - 1
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent + i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in start_count:
                                window_count += start_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in start_peaks:
                                    start_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                     extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['start_right'].add(identity)
    for position in sorted(ends, reverse=True):
        if position + upstream_buffer -1 not in end_peaks:
            window_count = 0
            for i in np.arange(0,10):
                window_position = position - i
                if window_position in end_count:
                    window_count += end_count[window_position]
            if window_count >= minimum_feature_count:
                original_bin=[]
                for shift in np.arange(-downstream_buffer, upstream_buffer):
                    end_peaks[position + shift] = position
                    original_bin.append(position + shift)
                if extend:
                    bins=[]
                    for i in np.arange(min(original_bin),max(original_bin)):
                        bin=0
                        for shift in np.arange(0,10):
                            if i+shift in end_count:
                                 bin+=end_count[i+shift]
                        bins.append(bin)

                    best_bin=max(bins)
                    extended = True
                    adjacent = position - downstream_buffer
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent - i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in end_count:
                                window_count += end_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in end_peaks:
                                    end_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent=adjacent_pos
                        if extended:
                            count_dict['end_left'].add(identity)

                    extended=True
                    adjacent=position + upstream_buffer - 1
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent + i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in end_count:
                                window_count += end_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in end_peaks:
                                    original_bin.append(element)
                                    end_peaks[element] = position
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['end_right'].add(identity)

    return start_peaks, end_peaks, count_dict


def sort_reads_into_splice_junctions(splice_dict, infile):
    start_end_dict, start_end_dict_mono = {}, {}
    for line in open(infile):
        a = line.strip().split('\t')
        read_chromosome, read_direction = a[13], a[8]
        name = a[9]
        show = False
        read_direction = '+'  # ignores read direction
        start, end = int(a[15]), int(a[16])
        if read_direction == '+':
            left_extra, right_extra = int(a[11]), int(a[10]) - int(a[12])
        if read_direction == '-':
            right_extra, left_extra = int(a[11]), int(a[10]) - int(a[12])

        failed = False
        identity = read_chromosome + '_'

        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]
        sequence=a[23]
        read=(name,sequence)
        for x in range(0, len(blocksizes) - 1):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            left_splice = blockstart + blocksize
            right_splice = int(blockstarts[x + 1])
            if right_splice - left_splice > 50:
                if read_chromosome not in splice_dict:
                    failed = True
                    break
                else:
                    if not splice_dict[read_chromosome].get(left_splice) or not splice_dict[read_chromosome].get(right_splice):
                        failed = True
                        break
                left_splice_site = splice_dict[read_chromosome][left_splice]
                right_splice_site = splice_dict[read_chromosome][right_splice]
                identity += str(left_splice_site) + '-' + str(right_splice_site) + '~'
        if not failed:
            if identity.split('_')[1] != '':
                if identity not in start_end_dict:
                    start_end_dict[identity] = []
                start_end_dict[identity].append((start, end,
                                                 read,
                                                 left_extra,
                                                 right_extra,
                                                 read_direction))
            else:
                if identity not in start_end_dict_mono:
                    start_end_dict_mono[identity] = []

                start_end_dict_mono[identity].append((start, end,
                                                      read,
                                                      left_extra,
                                                      right_extra,
                                                      read_direction))
    return start_end_dict, start_end_dict_mono


def group_mono_exon_transcripts(start_end_dict, start_end_dict_mono):
    for identity in start_end_dict_mono:
        previous_end = 0
        iso_counter = 0
        new_identity = identity + 'M' + str(iso_counter)
        positions = start_end_dict_mono[identity]
        for position in sorted(positions):
            start = position[0]
            end = position[1]
            if start > previous_end:
                iso_counter += 1
                new_identity = identity + 'M' + str(iso_counter)
                if new_identity not in start_end_dict:
                    start_end_dict[new_identity] = []
                start_end_dict[new_identity].append(position)
                previous_end = max(end, previous_end)
            else:
                if new_identity not in start_end_dict:
                    start_end_dict[new_identity] = []
                start_end_dict[new_identity].append(position)
                previous_end = end

    return start_end_dict


def define_start_end_sites(start_end_dict, start_end_dict_mono,upstream_buffer,downstream_buffer,minimum_feature_count):
    left_extras, right_extras = {}, {}
    isoform_counter, isoform_dict, IsoDict = 0, {}, {}

    start_end_dict = group_mono_exon_transcripts(start_end_dict, start_end_dict_mono)
    number_of_isoforms = len(start_end_dict)
    counter = 0
    count_dict={}
    count_dict['start_left']=set()
    count_dict['start_right']=set()
    count_dict['end_left']=set()
    count_dict['end_right']=set()
    count_dict['+']=set()
    count_dict['-']=set()

    for identity in sorted(start_end_dict):
        counter += 1
        starts, ends = [], []
        ends = []
        positions = start_end_dict[identity]
        length_of_positions = len(positions)
        indexes = np.random.choice(range(0, length_of_positions, 1),
                                   size=min(10000, length_of_positions),
                                   replace=False)
        subsampled_positions = []
        for index in indexes:
            subsampled_positions.append(positions[index])
        for position in subsampled_positions:
            starts.append(int(position[0]))
            ends.append(int(position[1]))

        start_dict, end_dict, count_dict = find_ends(starts, ends, identity, count_dict,upstream_buffer,downstream_buffer,minimum_feature_count)
        matched_positions = []
        left_extras[identity], right_extras[identity] = {}, {}
        for start, end, read, left_extra, right_extra, read_direction in positions:
            if int(start) in start_dict and int(end) in end_dict:
                left = start_dict[int(start)]
                right = end_dict[int(end)]
                if not left_extras[identity].get((left, right)):
                    left_extras[identity][(left, right)] = []
                    right_extras[identity][(left, right)] = []

                left_extras[identity][(left, right)].append(int(left_extra))
                right_extras[identity][(left, right)].append(int(right_extra))
                matched_positions.append((left, right, read, read_direction))

        median_left = {}
        median_right = {}
        for combination, values in left_extras[identity].items():
            median = np.median(values)
            median_left[combination] = median
        for combination, values in right_extras[identity].items():
            median = np.median(values)
            median_right[combination] = median

        for left, right, read, read_direction in matched_positions:
            medianLeft = median_left[(left, right)]
            medianRight = median_right[(left, right)]
            new_identity = identity + '_' + read_direction + '_' + str(left) \
                + '_' + str(right) + '_' \
                + str(round(medianLeft, 2)) \
                + '_' + str(round(medianRight, 2))
            if new_identity not in isoform_dict:
                isoform_counter += 1
                isoform_dict[new_identity] = isoform_counter

            IsoformName = str(isoform_dict[new_identity])
            if IsoformName not in IsoDict:
                IsoDict[IsoformName]=[]
            IsoDict[IsoformName].append(read)
    return IsoDict

def revComp(seq):
    '''Returns the reverse complement of a seq'''
    bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    return ''.join([bases[x] for x in list(seq)])[::-1]


def determine_consensus(reads,root,abpoa):
    '''Aligns and returns the consensus'''
    fasta_reads = []
    names=[]
    for read, seq in reads:
        fasta_reads.append((read, seq))
        names.append(read)

    length_fasta_reads=len(fasta_reads)

    indeces = np.random.choice(np.arange(0, length_fasta_reads),
                                   min(length_fasta_reads, 100), replace=False)

    indeces=set(indeces)
    subsample_fasta_reads = []
    lineCounter=-1
    for index in indeces:
        subsample_fasta_reads.append(fasta_reads[index])

    first = subsample_fasta_reads[0][1]
    sequences=[]
    seq_lengths=[]

    mm_align = mp.Aligner(seq=first, preset='map-ont')

    tmp_sequences=root+'.fasta'
    tmp_consensus=root+'.consensus.fasta'
    tmp_sequences_fh=open(tmp_sequences,'w')
    for read,sequence in subsample_fasta_reads:
        seq_lengths.append(len(sequence))
        for hit in mm_align.map(sequence):
             if hit.is_primary:
                 if hit.strand==-1:
                     sequence=mp.revcomp(sequence)
                 tmp_sequences_fh.write(f'>{read}\n{sequence}\n')
                 sequences.append(sequence)

    tmp_sequences_fh.close()

    if len(sequences)<=2:
         consensus_sequence = sequences[0]

    else:
        insert_length=np.median(seq_lengths)
        if insert_length<8000:
            os.system(f'{abpoa} -M 5 -r 0 {tmp_sequences} > {tmp_consensus} 2> apboa.messages')
        else:
            os.system(f'{abpoa} -M 5 -r 0 -S {tmp_sequences} > {tmp_consensus} 2> apboa.messages')

        consensus_sequence=''
        for consName,consSeq,consQ in mp.fastx_read(f'{tmp_consensus}'):
            consensus_sequence=consSeq
        if not consensus_sequence:
            consensus_sequence = sequences[0]

        os.system(f'rm {tmp_consensus}')
    os.system(f'rm {tmp_sequences}')

    res, mm_align, sequences = 0, 0, 0
    return consensus_sequence,names


#def split_consensus(sequences):

