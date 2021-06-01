#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import numpy as np
import argparse
import mappy as mp

def argParser():
    parser = argparse.ArgumentParser(
        description='Filters Isoforms generated by Mandalorion',
        add_help=True,
        prefix_chars='-',
    )

    parser.add_argument(
        '--path', '-p', type=str, action='store', default=os.getcwd(),
        help='Directory where all the files are/where they will end up.\
              Defaults to your current directory.',
    )
    parser.add_argument('--infile', '-i', type=str, action='store')
    parser.add_argument(
        '--config', '-c', type=str, action='store', default='',
        help='If you want to use a config file to specify paths to\
              programs, specify them here. Use for poa, racon, gonk,\
              blat, and minimap2 if they are not in your path.',
    )
    parser.add_argument('-n','--internal_ratio', type=float, action='store')
    parser.add_argument('-r', '--minimum_ratio', type=float)
    parser.add_argument('-R', '--minimum_reads', type=float)
    parser.add_argument('-G', '--genome_sequence', type=str)
    parser.add_argument('-a', '--adapter', type=str)
    parser.add_argument('-O', '--overhangs', type=str)
    parser.add_argument('-t', '--minimap2_threads', type=str)
    parser.add_argument('-A', '--Acutoff', type=str)
    parser.add_argument('-s', '--splice_window', type=str)
    parser.add_argument('-d', '--downstream_buffer', type=str)
    parser.add_argument('-I', '--minimum_isoform_length', type=str)
    parser.add_argument('-m', '--mandopath', type=str)
    parser.add_argument(
        '-e', '--ends', type=str, default='ATGGG,AAAAA',
        help='Ends of your sequences. Defaults to Smartseq ends.\
              Format: 5prime,3prime',
    )

    return vars(parser.parse_args())


def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, racon, consensus, blat, and emtrey
    possible = set(['minimap2', 'consensus', 'racon', 'blat', 'emtrey'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
    # check for missing programs
    # if missing, default to path
    for missing in possible - inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing) + ' from your path, not the config file.\n')
    return progs


args = argParser()
path = args['path']
infile = args['infile']
minimum_ratio = args['minimum_ratio']
minimum_reads = args['minimum_reads']
minimum_isoform_length = int(args['minimum_isoform_length'])
internal_ratio = args['internal_ratio']

genome = args['genome_sequence']
configIn = args['config']
adapter = args['adapter']
Acutoff = float(args['Acutoff'])
overhangs = np.array(args['overhangs'].split(','), dtype=int)
minimap2_threads = args['minimap2_threads']
ends = args['ends']
sw = int(args['splice_window'])
downstream_buffer = int(args['downstream_buffer'])
MandoPath=args['mandopath']

if args['config']:
    progs = configReader(args['config'])
    minimap2 = progs['minimap2']
    racon = progs['racon']
    consensus = progs['consensus']
    emtrey = progs['emtrey']
    blat = progs['blat']
else:
    minimap2, racon, emtrey, blat = 'minimap2', 'racon', 'emtrey', 'blat'
    consensus = 'consensus.py'

consensus = 'python3 ' + consensus

out2 = open(path + '/Isoform_Consensi_filtered.fasta', 'w')
out3 = open(path + '/Isoform_Consensi_filtered.aligned.out.clean.psl', 'w')
polyAWhiteListFile=path + '/polyAWhiteList.bed'

def reverse_complement(sequence):
    '''Returns the reverse complement of a sequence'''
    bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]


def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for name,seq,qual in mp.fastx_read(inFile):
        readDict[name]=seq
    return readDict


def get_count(isoform_list, chromosome, psl_dict):
    count = {}
    for isoform in isoform_list:
        info = psl_dict[isoform]
        coordinates = info[0]
        direction = info[1]
        number = int(isoform.split('_')[-1])
        start, end = coordinates[0], coordinates[-1]
        for base in np.arange(round(start, -1), round(end, -1), 10):
            if not count.get(chromosome + '_' + direction):
                count[chromosome + '_' + direction] = {}
            if not count[chromosome + '_' + direction].get(base):
                count[chromosome + '_' + direction][base] = number
            else:
                count[chromosome + '_' + direction][base] += number
    return count


def filter_isoforms(count, isoform_names, chromosome, psl_info, overhangs, minimum_isoform_length):
    filtered_isoforms = []
    for isoform in sorted(isoform_names):
        coverage_list = []
        info = psl_info[isoform]
        direction = info[8]
        start = int(info[15])
        end = int(info[16])

        number = int(isoform.split('_')[-1])
        for base in np.arange(round(start, -1), round(end, -1), 10):
            coverage_list.append(count[chromosome + '_' + direction][base])
        max_coverage = max(coverage_list)
        if (number / max_coverage) >= minimum_ratio:
            filtered_isoforms.append(isoform)
        else:
            sys.stderr.write(
                isoform + ' filtered because it at ' + str(number)
                + ' reads it only reached a ' + str(number / max_coverage)
                + ' ratio of expression in its locus which is below the minimum ratio of '
                + str(minimum_ratio) + '\n'
            )

    return filtered_isoforms


def look_for_contained_isoforms(isoform_list, chromosome, psl_dict, psl_info, genome_sequence,polyAWhiteList):
    internal_buffer = 20
    filtered_isoforms = []
    covered = {}
    covered[chromosome] = {}
    covered[chromosome]['+'] = {}
    covered[chromosome]['-'] = {}
    for isoform in isoform_list:
        info = psl_dict[isoform]
        counter = 0
        coordinates = info[0]
        direction = info[1]
        for coord in coordinates:
            counter += 1
            if counter % 2 == 1:
                start = coord
            if counter % 2 == 0:
                end = coord
                for base in range(start - sw, end + sw, 1):
                    if base not in covered[chromosome][direction]:
                        covered[chromosome][direction][base] = set()
                    covered[chromosome][direction][base].add(isoform)

    for isoform in isoform_list:
        info = psl_dict[isoform]
        status = set()
        coordinates = list(info[0])
        start, end = coordinates[0], coordinates[-1]
        coordinates[0] = min(coordinates[0] + internal_buffer, coordinates[1])
        coordinates[-1] = max(coordinates[-1] - internal_buffer, coordinates[-2])
        direction = info[1]
        if direction == '+':
            Acontent = genome_sequence[chromosome][end:end + 15].upper().count('A') / 15
            polyArange = np.arange(end + 3, end + 23, 1)
            polyApos=end
        elif direction == '-':
            Acontent = genome_sequence[chromosome][start - 15:start].upper().count('T') / 15
            polyArange = np.arange(start - 23, start - 3, 1)
            polyApos=start
        extend = set()
        extendDict = {}
        for base in polyArange:
            if base in covered[chromosome][direction]:
                for element in covered[chromosome][direction][base]:
                    if element not in extendDict:
                        extendDict[element] = 0
                    extendDict[element] += 1
        for element, value in extendDict.items():
            if value >= 10:
                extend.add(element)

        for isoform1 in psl_dict:
            status.add(isoform1)
        for coord in coordinates:
            counter += 1
            if counter % 2 == 1:
                start = coord
            if counter % 2 == 0:
                end = coord
                for base in range(start, end, 1):
                    status = status & covered[chromosome][direction][base]

        minimum_overlap = len(status) + len(extend)
        if minimum_overlap == 1:
            filtered_isoforms.append(isoform)
        else:
            show = False
            decision = True
            if '195088' in isoform:
                show = True
            if show:
                print('status', status)
                print('extend', extend)
            if len(extend) > 0:
                if Acontent > Acutoff:
                    if polyApos in polyAWhiteList[direction]:
                        sys.stderr.write(
                            isoform + ' would have been filtered because at least one isoform (including '
                            + str(list(extend)[0])
                            + ') is extending beyond its polyA site and the genomic A content at its putative polyA site is '
                            + str(Acontent) + ' which is higher than the cutoff set to ' + str(Acutoff)
                            + 'but it was kept because its polyA site was part of the polyA site whitelist\n'
                    )
                    else:
                        decision=False
                        sys.stderr.write(
                            isoform + ' filtered because at least one isoform (including '
                            + str(list(extend)[0])
                            + ') is extending beyond its polyA site and the genomic A content at its putative polyA site is '
                            + str(Acontent) + ' which is higher than the cutoff set to ' + str(Acutoff) + '\n'
                        )

            if decision == True:
                isoform_abundance = int(isoform.split('_')[-1])
                for match in status:
                    if show:
                        print(match)

                    if match != isoform:
                        match_abundance = int(match.split('_')[-1])
                        match_coordinates = psl_dict[match][0]
                        duplicate_dict = {}
                        for junction in range(1, len(match_coordinates) - 1, 2):
                            junc1 = match_coordinates[junction]
                            junc2 = match_coordinates[junction + 1]
                            for base1 in range(junc1 - sw, junc1 + sw, 1):
                                duplicate_dict[base1] = {}
                                for base2 in range(junc2 - sw, junc2 + sw, 1):
                                    duplicate_dict[base1][base2] = 1

                        matches = set()
                        all_junctions = set()
                        for junction in range(1, len(coordinates) - 1, 2):
                            all_junctions.add(junction)
                            junc1 = coordinates[junction]
                            junc2 = coordinates[junction + 1]
                            for base1 in range(junc1 - sw, junc1 + sw, 1):
                                if base1 in duplicate_dict:
                                    for base2 in range(junc2 - sw, junc2 + sw, 1):
                                        if base2 in duplicate_dict[base1]:
                                            matches.add(junction)
                        if show:
                            print('pre', all_junctions, matches)
                            print(coordinates, match_coordinates)
                        if len(all_junctions - matches) > 0:
                            if show:
                                print('check', matches, all_junctions)
                            continue

                        elif len(all_junctions - matches) == 0:
                            coordinates = list(info[0])
                            if show:
                                print((isoform_abundance / match_abundance), internal_ratio)
                            if (isoform_abundance / match_abundance) < internal_ratio:
                                sys.stderr.write(
                                    isoform
                                    + ' filtered because it is internal to (all bases and splice junctions contained in) '
                                    + match + ' and expressed at ' + str(isoform_abundance) + ' reads compared to '
                                    + str(match_abundance)
                                    + ' reads for the isoform containing it which is below that internal ratio of '
                                    + str(internal_ratio) + '\n'
                                )
                                decision=False
                                break

                            elif np.abs(coordinates[0] - match_coordinates[0]) < downstream_buffer:
                                if np.abs(coordinates[-1] - match_coordinates[-1]) < downstream_buffer:
                                    if isoform_abundance < match_abundance:
                                        sys.stderr.write(isoform + ' filtered because it is internal (all bases and splice junctions contained in) and almost identical to ' + match + '\n')
                                        decision=False
                                        break

            if show:
                print(decision)
            if  decision==True:
                filtered_isoforms.append(isoform)

    return filtered_isoforms


def read_splice_file(SS_file):
    splice_dict = {}
    for line in open(SS_file):
        a = line.strip().split('\t')
        start = int(a[1])
        end = int(a[2])
        name = a[3].split('_')[0]
        splice_dict[name] = (start, end)
    return splice_dict


def simplify(infile, outfile, namefile):
    isoforms = read_fasta(infile)
    out = open(outfile, 'w')
    out1 = open(namefile, 'w')
    counter = 0
    for isoform, sequence in isoforms.items():
        if len(sequence) > 0:
            counter += 1
            out.write('>Isoform' + '_' + str(counter) + '_' + isoform.split('_')[-1] + '\n' + sequence + '\n')
            out1.write(isoform + '\tIsoform' + '_' + str(counter) + '\n')
    out.close()
    out1.close()


def parse_clean_psl(psl_file, target_chromosome):
    psl_dict = {}
    psl_info = {}
    isoform_list = []
    first_alignment = {}
    for line in open(psl_file):
        a = line.strip().split('\t')

        chromosome = a[13]
        if chromosome != target_chromosome:
            continue
        readstart = int(a[11])
        readend = int(a[12])
        readlength = readend - readstart
        direction = a[8]

        if direction == '+':
            overhang5 = int(a[11])
            overhang3 = int(a[10]) - int(a[12])
        if direction == '-':
            overhang3 = int(a[11])
            overhang5 = int(a[10]) - int(a[12])

        name = a[9]
        abundance = int(a[9].split('_')[-1])
        if readlength >= minimum_isoform_length:
            if abundance >= minimum_reads:
                if overhangs[0] <= overhang5 <= overhangs[1] and overhangs[2] <= overhang3 <= overhangs[3]:
                    if name not in first_alignment:
                        isoform_list.append(name)
                        psl_info[name] = a
                        blocksizes = a[18].split(',')[:-1]
                        blockstarts = a[20].split(',')[:-1]
                        psl_dict[name] = [[], direction]
                        for index in np.arange(0, len(blocksizes), 1):
                            blockstart = int(blockstarts[index])
                            blocksize = int(blocksizes[index])
                            blockend = blockstart + blocksize
                            psl_dict[name][0].append(blockstart)
                            psl_dict[name][0].append(blockend)
                else:
                    sys.stderr.write(
                        name + ' filtered because at ' + str(overhang5) + ' and ' + str(overhang3)
                        + ' its number of overhanging bases did not fall within the predefined bins of '
                        + str(overhangs[0]) + '-' + str(overhangs[1]) + ' and ' + str(overhangs[2])
                        + '-' + str(overhangs[3]) + '\n'
                    )

            else:
                sys.stderr.write(
                    name + ' filtered because it at ' + str(abundance)
                    + ' reads it did not match the minimum absolute read requirement of '
                    + str(minimum_reads) + '\n'
                )
        else:
            sys.stderr.write(
                name + ' filtered because at ' + str(readlength)
                + 'nt it did not match the minimum isoform length requirement of '
                + str(minimum_isoform_length) + '\n'
            )
        first_alignment[name] = 1

    return psl_dict, psl_info, set(isoform_list)


def collect_chromosomes(isoform_psl):
    chromosomes = set()
    for line in open(isoform_psl):
        a = line.strip().split('\t')
        chromosome = a[13]
        chromosomes.add(chromosome)
    chromosomes = sorted(list(chromosomes))
    return chromosomes


def write_isoforms(isoform_list, isoforms, psl_info):
    for isoform in isoform_list:
        info = psl_info[isoform]
        sequence = isoforms[isoform]
        out2.write('>%s\n%s\n' % (isoform, sequence))
        out3.write('\t'.join(info) + '\n')

def readWhiteList(polyA,chromosome):

    WhiteList={}
    WhiteList['+']=set()
    WhiteList['-']=set()

    for line in open(polyA):
        a=line.strip().split('\t')
        if chromosome == a[0]:
            for pos in np.arange(int(a[1]),int(a[2]),1):
                WhiteList[a[5]].add(pos)
    return WhiteList


def main(infile):
    print('simplifying isoform names')
    temp_fasta = path + '/isoform_tmp.fasta'
    simplify(infile, temp_fasta, path + '/Isoform_long_names.txt')
    genome_sequence = read_fasta(genome)
    os.system('python3 %s/%s -i %s -a %s -o %s -c %s -e %s' % (MandoPath,'postprocessingIsoforms.py', temp_fasta, adapter, path, configIn, ends))
    print('reading fasta')
    processed_isoforms = path + 'Isoforms_full_length_consensus_reads.fasta'
    isoforms = read_fasta(processed_isoforms)
    sam_file = path + '/Isoforms.aligned.out.sam'
    psl_file = path + '/Isoforms.aligned.out.psl'
    clean_psl_file = path + '/Isoforms.aligned.out.clean.psl'
    print('aligning reads')
    print('%s -G 400k --secondary=no -ax splice:hq -t %s %s %s > %s ' % (minimap2, minimap2_threads, genome, processed_isoforms, sam_file))
    os.system('%s -G 400k --secondary=no -ax splice:hq -t %s %s %s > %s ' % (minimap2, minimap2_threads, genome, processed_isoforms, sam_file))
    os.system('%s -i %s > %s ' % (emtrey, sam_file, psl_file))
    os.system('python3 %s/%s %s %s ' % (MandoPath,'clean_psl.py', psl_file, clean_psl_file))
    print('collecting chromosomes')
    chromosomes = collect_chromosomes(clean_psl_file)
    for chromosome in chromosomes:
        print(chromosome)
        sys.stderr.write(chromosome + '\n')
        print('reading polyA white list')
        polyAWhiteList=readWhiteList(polyAWhiteListFile,chromosome)
        print('reading in isoforms and applying absolute filters for abundance, lengths, and overhangs')
        psl_dict, psl_info, isoform_list = parse_clean_psl(clean_psl_file, chromosome)
        print('getting isoform loci read counts')
        count = get_count(isoform_list, chromosome, psl_dict)
        print('filtering isoforms for relative read coverage starting with', len(isoform_list), 'isoforms')
        isoform_list = filter_isoforms(count, isoform_list, chromosome, psl_info, overhangs, minimum_isoform_length)
        print('finding fully contained isoforms in', len(isoform_list), 'remaining isoforms')
        isoform_list = look_for_contained_isoforms(isoform_list, chromosome, psl_dict, psl_info, genome_sequence,polyAWhiteList)
        print('writing', len(isoform_list), 'isoforms to file')
        write_isoforms(isoform_list, isoforms, psl_info)


main(infile)
