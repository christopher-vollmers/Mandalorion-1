#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import argparse
import mappy as mp
import time
from time import localtime, strftime

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/utils/'
sys.path.append(os.path.abspath(PATH))

import SpliceDefineConsensus


VERSION = "v4.5.0 - Somehow Isoform Returned"

parser = argparse.ArgumentParser(usage='\n\nRunning with default parameters:\n\npython3 Mando.py -p . -g gencodeV29.gtf -G hg38.fasta -f Consensus_reads_noAdapters_noPolyA_5->3.fofn\n')

parser.add_argument(
    '-p', '--path', type=str, help='Directory to put output files into'
)
parser.add_argument(
    '-u', '--upstream_buffer', type=str, default='10',
    help='''Defines upstream leniency window for
            polyA and TSS determination (default 10)'''
)
parser.add_argument(
    '-d', '--downstream_buffer', type=str, default='50',
    help='''Defines downstream leniency window for polyA
            and TSS determination (default 50)'''
)
parser.add_argument(
    '-g', '--genome_annotation', type=str, default='None',
    help='''Genome annotation file (gtf).
            Is used to identify individual annotated splice sites.
            If -W is set it is also used to white-list annotated polyA sites'''
)
parser.add_argument(
    '-G', '--genome_sequence', type=str, help='Genome file (fasta)'
)
parser.add_argument(
    '-r', '--minimum_ratio', type=str, default='0.01',
    help='''Proportion of reads that align to a locus
            required for an isoform (default 0.01)'''
)
parser.add_argument('-i', '--minimum_internal_ratio', type=str, default='1')
parser.add_argument(
    '-R', '--minimum_reads', type=str, default='3',
    help='''Minimum number of reads to make an isoform (default 3)'''
)

parser.add_argument(
    '-f', '--Consensus_reads', type=str,
    help='''Fasta/fastq file with R2C2/PacBio consensus reads,
            can be entered as a single file path,
            a comma separated list of  file paths,
            or a path to a file of filenames file (has to end on .fofn) that contains one file path per line'''
)

parser.add_argument(
    '-O', '--overhangs', type=str, default='0,40,0,40',
    help='''Defines bounds for unaligned bases on ends. Format:
            min5prime,max5prime,min3prime,max3prime (default 0,40,0,40)'''
)
parser.add_argument(
    '-t', '--minimap2_threads', type=str, default='8',
    help='Number of threads to use when running minimap and consensus calling (default 8)'
)

parser.add_argument(
    '-I', '--minimum_isoform_length', type=str, default='200',
    help='Minimum length in nt of isoforms that will be considered (default 200)'
)
parser.add_argument(
    '-n', '--minimum_feature_count', type=str, default='2',
    help='''Features (starts,ends, new splice sites) will be considered
            if they are present in this number of reads (default 2)'''
)
parser.add_argument(
    '-w', '--splice_site_window', type=str, default='1',
    help='''reads spliced within this number of nucleotides on each side
            of a splice site will be considered spliced at this site (default 1)'''
)
parser.add_argument(
    '-A', '--Acutoff', type=str, default='0.5',
    help='''Isoforms with A content of more than this cutoff in a 30nt
            window surrounding their polyA site will be discarded (default 0.5)'''
)
parser.add_argument(
    '-W', '--white_list_polyA', type=str, default='0',
    help='''If set, polyA sites that fall within +/-20nt of annotated transcript ends will not be filtered regardless of Acutoff set with -A.
            Annotated transcript ends will be  taken from annotation file given with -g.
            Only transcripts having one of the provided comma separated values in the line of their exon features will be used.
            E.g. setting -W to [SIRV,"basic"] will whitelist spike-in SIRV transcripts and transcripts considered "basic", i.e. high confidence full length, in the gencode annotation.
            Setting -W to [exon] should include all transcripts in the gtf file.
            This is feature only checks whether the provided values are in the line, not where they are.
            That means that setting -W to [chr1] will whitelist all transcripts on chromosome 1 but also transcripts with "chr1" in their name, so it's a bit dangerous'''
)


parser.add_argument(
    '-m', '--multi_exon_only', default='0',action='store_const', const='1',
    help='''If used, Mandalorion will filter all single exon isoforms'''
)

parser.add_argument(
    '-j', '--junctions', default='gtag,gcag,atac,ctac,ctgc,gtat',action='store',type=str,
    help='''base context required (in lower-case for either strand) to call unannotated splice junctions (default gtag,gcag,atac,ctac,ctgc,gtat)'''
)

parser.add_argument(
    '-M', '--Modules', default='APDFQ',
    help='''Defines what modules of Mandalorion will be run. By default this includes:
            A - Alignment,
            P - .sam to .clean.psl conversion,
            D - defining isoforms,
            F - Filtering isoforms,
            Q - Quantifying isoforms.
            Each module needs the output of the modules run before it to function properly.
            Running individual modules can be useful if you want to for example filter with different parameters without rerunning the whole pipeline. (default APDFQ)'''
)


parser.add_argument(
    '-v', '--version', action='version', version=VERSION,
    help='Prints Mandalorion version'
)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()

path = args.path + '/'  # path where you want your output files to go
temp_path = path + '/tmp/'
upstream_buffer = args.upstream_buffer
downstream_buffer = args.downstream_buffer
genome_annotation = args.genome_annotation
genome_sequence = args.genome_sequence
minimum_ratio = args.minimum_ratio
minimum_internal_ratio = args.minimum_internal_ratio
minimum_reads = args.minimum_reads
fasta_files = args.Consensus_reads
overhangs = args.overhangs
minimap2_threads = args.minimap2_threads
minimum_isoform_length = args.minimum_isoform_length
window = args.splice_site_window
feature_count = args.minimum_feature_count
Acutoff = args.Acutoff
white_list_polyA=args.white_list_polyA
multi_exon_only=args.multi_exon_only
junctions=args.junctions
Modules=args.Modules

MandoPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'

if '.fofn' in fasta_files:
    fastaList=[]
    for line in open(fasta_files):
        fasta=line.strip()
        fastaList.append(fasta)
else:
    fastaList=fasta_files.split(',')

if not os.path.isdir(path):
    os.system('mkdir %s' % path)
    command=open(path+'/Mando.log','w')
else:
    command=open(path+'/Mando.log','a')

command.write('\nMandalorion "'+ VERSION +'" was run on '\
               + strftime("%Y-%m-%d %H:%M:%S", localtime())+'\n'\
               + 'with the following parameters\n'\
               + str(args).replace('Namespace(','').replace(')','') + '\n')
command.close()
if not os.path.isdir(temp_path):
    os.system('mkdir %s' % temp_path)

minimap2 = MandoPath+'minimap2/minimap2'
emtrey = MandoPath+'/emtrey.py'
abpoa=MandoPath+'/abPOA-v1.4.1/bin/abpoa'

print('\n----------------------------------------------------------------\
       \nMandalorion - Isoform identification and quantification\
       \nVersion -', VERSION,'\
       \n----------------------------------------------------------------\n')


sam_file = temp_path + '/mm2Alignments.sam'
if 'A' in Modules:
    print('\n----------------------------\
           \n    Module A - Alignment\
           \n----------------------------\n')
    tempFasta=temp_path+'/Combined.fasta'

    print('\tchecking input fasta/q files')
    reads=False
    if len(fastaList)==1:
        print('\t1 fasta/q file provided')
        use_file=False
        fasta=fastaList[0]
        if os.path.exists(fasta) and os.path.getsize(fasta)>0:
            use_file=True
        if use_file:
            tempFasta=fastaList[0]
            reads=True
        else:
            print('\t',fasta, 'does not exist or is an empty file')
    else:
        print('\tcombining '+str(len(fastaList))+' input fasta/q files')
        out=open(tempFasta,'w')
        for fasta in fastaList:
            use_file=False
            if os.path.exists(fasta) and os.path.getsize(fasta)>0:
                use_file=True
            if use_file:
                reads=True
                for name,seq,qual in mp.fastx_read(fasta):
                    out.write('>%s\n%s\n' %(name,seq))
            else:
                print('\t',fasta, 'does not exist or is an empty file')
        out.close()
    if reads:
        os.system(
            '%s -G 400k --secondary=no -ax splice:hq --cs=long -uf -t %s %s %s > %s '
            % (minimap2, minimap2_threads, genome_sequence, tempFasta, sam_file)
        )
    else:
        print('\t no reads were provided. Alignment will not be performed')

psl_file = temp_path + '/mm2Alignments.psl'
clean_psl_file = temp_path + '/mm2Alignments.clean.psl'
clean_sorted_psl_file = temp_path + '/mm2Alignments.clean.sorted.psl'
if 'P' in Modules:
    print('\n------------------------------------------------------\
           \n    Module P - .sam to .clean.sorted.psl conversion\
           \n------------------------------------------------------\n')

    input=False
    if os.path.exists(sam_file) and os.path.getsize(sam_file)>0:
        input=True

    if input:
        print('\tconverting sam output to psl format')
        os.system('python3 %s -i %s -o %s -m -t %s' % (emtrey, sam_file, psl_file,minimap2_threads))
        print('\tcleaning psl file of small Indels')
        SpliceDefineConsensus.clean_psl(psl_file, clean_psl_file,True)
        print('\tsorting clean psl file')
        os.system('sort -T %s -k 14,14 -k 16,17n %s > %s' %(temp_path,clean_psl_file,clean_sorted_psl_file))
        print('\treading and splitting psl file into loci')
        if os.path.isdir(f'{temp_path}tmp_SS/'):
            os.system('rm -r %s/%s' % (temp_path,'tmp_SS/'))
        os.system('mkdir %s/%s' % (temp_path,'tmp_SS/'))
        SpliceDefineConsensus.get_loci(clean_sorted_psl_file,temp_path+'/tmp_SS',minimum_reads)
    else:
        print('\tno or empty SAM file was provided. File conversions and parsing not performed')


if 'D' in Modules:
    print('\n----------------------------------------\
           \n      Module D - defining isoforms\
           \n----------------------------------------\n')
    everything_present=True

    if not os.path.exists(clean_sorted_psl_file) or os.path.getsize(clean_sorted_psl_file)==0:
        print('\tclean sorted psl file missing or empty')
        everything_present=False
    for fasta_file in fastaList:
        if not os.path.exists(fasta_file) or os.path.getsize(fasta_file)==0:
            print('\t',fasta_file,'missing or empty')
            everything_present=False

    if everything_present:
        os.system(
            'python3 %s/defineIsoforms.py -i %s -p %s -c %s -g %s -w %s -m %s -W %s -n %s -j %s -u %s -d %s -a %s'
            % (
                MandoPath,
                clean_sorted_psl_file,
                temp_path,
                '0.1',
                genome_annotation,
                window,
                feature_count,
                white_list_polyA,
                minimap2_threads,
                junctions,
                upstream_buffer,
                downstream_buffer,
                abpoa
            )
        )
    else:
        print('\tone or more input files missing or empty. Isoforms not defined')

if 'F' in Modules:
    print('\n-------------------------------------\
           \n    Module F - filtering isoforms\
           \n-------------------------------------\n')

    everything_present=True
    if not os.path.exists(temp_path + '/Isoform_Consensi.fasta') or os.path.getsize(temp_path + '/Isoform_Consensi.fasta')==0:
        print('\tisoforms fasta missing or empty')
        everything_present=False

    if not os.path.exists(genome_sequence) or os.path.getsize(genome_sequence)==0:
        print('\tgenome sequence fasta missing or empty')
        everything_present=False

    if everything_present:
        os.system(
            'python3 %s/filterIsoforms.py \
                -p %s -i %s -r %s -R %s -n %s -G %s \
                -O %s -t %s -A %s -s %s -d %s -I %s -m %s -M %s 2> %s'
            % (
                MandoPath,
                temp_path,
                temp_path + '/Isoform_Consensi.fasta',
                minimum_ratio,
                minimum_reads,
                minimum_internal_ratio,
                genome_sequence,
                overhangs,
                minimap2_threads,
                Acutoff,
                window,
                downstream_buffer,
                minimum_isoform_length,
                MandoPath,
                multi_exon_only,
                temp_path + '/filter_reasons.txt',
            )
        )

        os.system('sort -T %s -k 9,9 -k 14,14 -k 16,17n %s > %s' %(temp_path, temp_path + '/Isoforms.filtered.clean.psl',temp_path + '/Isoforms.sorted.psl'))
        os.system('mv %s %s' %(temp_path + '/Isoforms.sorted.psl',temp_path + '/Isoforms.filtered.clean.psl'))

        print('\tgrouping isoforms, assigning them to genes (if annotation is provided), and creating gtf')
        os.system(
            'python3 %s/groupIsoforms.py -i %s -o %s -g %s -t %s'
                  % (MandoPath,
                     temp_path + '/Isoforms.filtered.clean.psl',
                     temp_path + '/Isoforms.filtered.clean.genes',
                     genome_annotation,
                     temp_path + '/Isoforms.filtered.clean.gtf'
            )
        )
        os.system('scp ' + temp_path + '/Isoforms.filtered.* ' + path)
    else:
        print('\tone or more input files missing or empty. Isoforms not filtered')

if 'Q' in Modules:
    print('\n---------------------------------------\
           \n    Module Q - quantifying isoforms\
           \n---------------------------------------\n')
    os.system(
        'python3 %s/assignReadsToIsoforms.py -m %s -f %s'
        % (MandoPath,temp_path, fasta_files)
    )
    os.system('scp ' + temp_path + '/*.quant ' + path)
    os.system('scp ' + temp_path + '/*.rpm ' + path)
    os.system('scp ' + temp_path + '/*.frac ' + path)
