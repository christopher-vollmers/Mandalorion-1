import sys
import argparse
import mappy as mp
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mandalorion_output_folder', type=str)
parser.add_argument('-f', '--fasta_files', type=str, help='comma separate list of fasta file locations')



args = parser.parse_args()
mandalorion_folder=args.mandalorion_output_folder
fasta_files=args.fasta_files
filtered_isoforms=mandalorion_folder+'/Isoforms.filtered.clean.psl'
gene_file=mandalorion_folder+'/Isoforms.filtered.clean.genes'
r2i=mandalorion_folder+'/reads2isoforms.txt'


outq = open(mandalorion_folder+'/Isoforms.filtered.clean.quant','w')
outtpm = open(mandalorion_folder+'/Isoforms.filtered.clean.rpm','w')
outfrac = open(mandalorion_folder+'/Isoforms.filtered.clean.frac','w')
outStartQ = open(mandalorion_folder+'/TSSs.quant','w')
outStartTpm = open(mandalorion_folder+'/TSSs.rpm','w')
outStartFrac = open(mandalorion_folder+'/TSSs.frac','w')
outEndQ = open(mandalorion_folder+'/PolyAs.quant','w')
outEndTpm = open(mandalorion_folder+'/PolyAs.rpm','w')
outEndFrac = open(mandalorion_folder+'/PolyAs.frac','w')
outJunctionQ = open(mandalorion_folder+'/Junctions.quant','w')
outJunctionTpm = open(mandalorion_folder+'/Junctions.rpm','w')
outJunctionFrac = open(mandalorion_folder+'/Junctions.frac','w')

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for name,seq,qual in mp.fastx_read(inFile):
         readDict[name] = ''

    return readDict

def find_regions(positions):

    locDict={}
    for chr_dir,posList in positions.items():
        sortedPosList=sorted(posList)
        isoforms=set()
        locs=[0]
        for pos,isoform in sortedPosList:
            if pos>max(locs)+100:
                dist2=pos-max(locs)
                region=f'{chr_dir}~{min(locs)}~{max(locs)+min(100,int(dist2/2)-1)}'
                for isoform1 in isoforms:
                    locDict[isoform1]=region
                locs=[pos-min(100,int(dist2/2)-1),pos]
                isoforms=set()
                isoforms.add(isoform)
            else:
                isoforms.add(isoform)
                locs.append(pos)
        if locs:
            region=f'{chr_dir}~{min(locs)}~{max(locs)+100}'
            for isoform1 in isoforms:
                locDict[isoform1]=region
    return locDict

def get_features(filtered_isoforms):
    starts={}
    ends={}
    isoJunctions={}
    for line in open(filtered_isoforms):
        a=line.strip().split('\t')
        isoform,chromosome,direction=a[9],a[13],a[8]
        if direction=='-':
            start,end=int(a[16]),int(a[15])
        else:
            start,end=int(a[15]),int(a[16])
        if chromosome+'~'+direction not in starts:
            starts[chromosome+'~'+direction]=[]
            ends[chromosome+'~'+direction]=[]
        starts[chromosome+'~'+direction].append((start,isoform))
        ends[chromosome+'~'+direction].append((end,isoform))
        blocksizes=np.array(a[18].split(',')[:-1],dtype=int)
        blockstarts=np.array(a[20].split(',')[:-1],dtype=int)
        isoJunctions[isoform]=[]
        if len(blockstarts)>1:
            for i in range(0,len(blockstarts)-1,1):
                junction_start=blockstarts[i]+blocksizes[i]
                junction_end=blockstarts[i+1]
                isoJunctions[isoform].append(f'{chromosome}~{direction}~{junction_start}~{junction_end}')

    isoStarts=find_regions(starts)
    isoEnds=find_regions(ends)

    return isoJunctions,isoStarts,isoEnds

def read_filtered_isoforms(filtered_isoforms,r2i_dict,sampleList,readMapDict,isoformReadCounts,totalReadCounts,geneDict):
    isoJunctions,isoStarts,isoEnds = get_features(filtered_isoforms)
    juncDict={}
    startDict={}
    endDict={}

    geneQuant={}
    for line in open(filtered_isoforms):
        a=line.strip().split('\t')
        isoform=a[9]
        gene=geneDict[isoform]
        if gene not in geneQuant:
            geneQuant[gene]={}
            for sample in sampleList:
                geneQuant[gene][sample]=0
        for name in r2i_dict[isoform]:
            sample=readMapDict[name]
            geneQuant[gene][sample]+=1

    start2gene={}
    end2gene={}
    junction2gene={}
    for line in open(filtered_isoforms):
        a=line.strip().split('\t')
        isoform=a[9]
        gene=geneDict[isoform]
        junctions = isoJunctions[isoform]
        start=isoStarts[isoform]
        end=isoEnds[isoform]
        start2gene[start]=gene
        end2gene[end]=gene
        quantDict={}


        for sample in sampleList:
            quantDict[sample]=0

        for junction in junctions:
            junction2gene[junction]=gene
            if junction not in juncDict:
                juncDict[junction]={}
                for sample in sampleList:
                    juncDict[junction][sample]=0

        if start not in startDict:
            startDict[start]={}
            for sample in sampleList:
                startDict[start][sample]=0

        if end not in endDict:
            endDict[end]={}
            for sample in sampleList:
                endDict[end][sample]=0

        for name in r2i_dict[isoform]:
            sample=readMapDict[name]
            quantDict[sample]+=1

            endDict[end][sample]+=1
            startDict[start][sample]+=1
            for junction in junctions:
                juncDict[junction][sample]+=1


        outq.write(a[9]+'\t'+gene+'\t')
        outtpm.write(a[9]+'\t'+gene+'\t')
        outfrac.write(a[9]+'\t'+gene+'\t')
        for sample in sampleList:
            value=quantDict[sample]
            isoformReads=isoformReadCounts[sample]
            geneReads=geneQuant[gene][sample]
            totalReads=totalReadCounts[sample]
            outq.write(str(value)+'\t')
            outtpm.write(str(round((value/totalReads)*1000000,3))+'\t')
            outfrac.write(f'{value}/{geneReads}\t')
        outq.write('\n')
        outtpm.write('\n')
        outfrac.write('\n')
    keys=list(startDict.keys())
    sortedKeys=sorted(keys, key= lambda x: (x.split('~')[1],x.split('~')[0],int(x.split('~')[2])))

    for start in sortedKeys:
        gene=start2gene[start]
        Chromosome,Direction,Start,End=start.split('~')
        outStartQ.write(f'{Chromosome}\t{Start}\t{End}\t{start}\t1000\t{Direction}\t{gene}\t')
        outStartTpm.write(f'{Chromosome}\t{Start}\t{End}\t{start}\t1000\t{Direction}\t{gene}\t')
        outStartFrac.write(f'{Chromosome}\t{Start}\t{End}\t{start}\t1000\t{Direction}\t{gene}\t')
        for sample in sampleList:
            value=startDict[start][sample]
            totalReads=totalReadCounts[sample]
            geneReads=geneQuant[gene][sample]
            outStartQ.write(str(value)+'\t')
            outStartTpm.write(str(round((value/totalReads)*1000000,3))+'\t')
            outStartFrac.write(f'{value}/{geneReads}\t')
        outStartQ.write('\n')
        outStartTpm.write('\n')
        outStartFrac.write('\n')

    keys=list(juncDict.keys())
    sortedKeys=sorted(keys, key= lambda x: (x.split('~')[1],x.split('~')[0],int(x.split('~')[2])))

    for junction in sortedKeys:
        gene=junction2gene[junction]
        Chromosome,Direction,Start,End=junction.split('~')
        outJunctionQ.write(f'{Chromosome}\t{Start}\t{End}\t{junction}\t1000\t{Direction}\t{gene}\t')
        outJunctionTpm.write(f'{Chromosome}\t{Start}\t{End}\t{junction}\t1000\t{Direction}\t{gene}\t')
        outJunctionFrac.write(f'{Chromosome}\t{Start}\t{End}\t{junction}\t1000\t{Direction}\t{gene}\t')
        for sample in sampleList:
            value=juncDict[junction][sample]
            totalReads=totalReadCounts[sample]
            geneReads=geneQuant[gene][sample]
            outJunctionQ.write(str(value)+'\t')
            outJunctionTpm.write(str(round((value/totalReads)*1000000,3))+'\t')
            outJunctionFrac.write(f'{value}/{geneReads}\t')

        outJunctionQ.write('\n')
        outJunctionTpm.write('\n')
        outJunctionFrac.write('\n')

    keys=list(endDict.keys())
    sortedKeys=sorted(keys, key= lambda x: (x.split('~')[1],x.split('~')[0],int(x.split('~')[2])))

    for end in sortedKeys:
        gene=end2gene[end]
        Chromosome,Direction,Start,End=end.split('~')
        outEndQ.write(f'{Chromosome}\t{Start}\t{End}\t{end}\t1000\t{Direction}\t{gene}\t')
        outEndTpm.write(f'{Chromosome}\t{Start}\t{End}\t{end}\t1000\t{Direction}\t{gene}\t')
        outEndFrac.write(f'{Chromosome}\t{Start}\t{End}\t{end}\t1000\t{Direction}\t{gene}\t')
        for sample in sampleList:
            value=endDict[end][sample]
            totalReads=totalReadCounts[sample]
            geneReads=geneQuant[gene][sample]
            outEndQ.write(str(value)+'\t')
            outEndTpm.write(str(round((value/totalReads)*1000000,3))+'\t')
            outEndFrac.write(f'{value}/{geneReads}\t')

        outEndQ.write('\n')
        outEndTpm.write('\n')
        outEndFrac.write('\n')



def mapReadLocation(fastaList):
    sampleList=[]
    readMapDict={}
    totalReadCounts={}
    for line in fastaList:
        location=line.strip()
        totalReadCounts[location]=0
        sampleList.append(location)
        reads=read_fasta(location)
        for name,seq,qual in mp.fastx_read(location):
            readMapDict[name]=location
            totalReadCounts[location]+=1

    for outFile in (outq,outtpm,outfrac):
        outFile.write('Isoform\tGene\t')
    for outFile in (outStartQ,outStartTpm,outStartFrac,outEndQ,outEndTpm,outEndFrac,outJunctionQ,outJunctionTpm,outJunctionFrac):
       outFile.write('Chr\tStart\tEnd\tName\tScore\tDirection\tGene\t')


    for sample in sampleList:
        for outFile in (outq,outtpm,outfrac,outStartQ,outStartTpm,outStartFrac,outEndQ,outEndTpm,outEndFrac,outJunctionQ,outJunctionTpm,outJunctionFrac):
            outFile.write(sample+'\t')
    for outFile in (outq,outtpm,outfrac,outStartQ,outStartTpm,outStartFrac,outEndQ,outEndTpm,outEndFrac,outJunctionQ,outJunctionTpm,outJunctionFrac):
        outFile.write('\n')

    return sampleList,readMapDict,totalReadCounts


def read_r2i(r2i,readMapDict):
    r2i_dict={}
    isoformReadCounts={}
    for line in open(r2i):
        a=line.strip().split('\t')
        read=a[0]
        isoform=a[1]
        location=readMapDict[read]
        if location not in isoformReadCounts:
            isoformReadCounts[location]=0
        if isoform not in r2i_dict:
            r2i_dict[isoform]=[]
        r2i_dict[isoform].append(read)
        isoformReadCounts[location]+=1
    return r2i_dict,isoformReadCounts

def getGenes(gene_file):
    genes=set()
    geneDict={}
    with open(gene_file) as f:
        for line in f:
            a=line.strip().split('\t')
            isoform=a[0]
            if len(a)>5:
                gene=a[5]
            else:
                gene=a[1]
            geneDict[isoform]=gene
    return geneDict,genes

if '.fofn' in fasta_files:
    fastaList=[]
    for line in open(fasta_files):
        fasta=line.strip()
        fastaList.append(fasta)
else:
    fastaList=fasta_files.split(',')




print('\tassigning isoforms to genes')
geneDict,genes=getGenes(gene_file)
print('\tmap reads to samples')
sampleList,readMapDict,totalReadCounts=mapReadLocation(fastaList)
print('\tmap reads to isoforms')
r2i_dict,isoformReadCounts=read_r2i(r2i,readMapDict)
print('\tquantify isoforms and isoform features')
read_filtered_isoforms(filtered_isoforms,r2i_dict,sampleList,readMapDict,isoformReadCounts,totalReadCounts,geneDict)
