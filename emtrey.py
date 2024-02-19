import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--threads','-t',default='50',type=str,action='store',help='number of threads to be used')
parser.add_argument('--outputFile','-o',type=str,action='store',help='output file goes here')
parser.add_argument('--inputFile','-i',type=str,action='store',help='path to input file')
parser.add_argument('--batch','-b',default='100000',type=str,action='store',help='lines are processed in batches of this size')
parser.add_argument('--mando','-m',default=False,action='store_true',help='generates invalid psl file for internal Mandalorion use')



args = parser.parse_args()

threads=int(args.threads)
inFile=args.inputFile
outputFile=args.outputFile
mandoBool=args.mando
batch=int(args.batch)

if mandoBool:
    mando='-m'
else:
    mando=''

MandoPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'

out=open(outputFile,'w')
out.close()

mm=True



def delegatingPslConversion(samFile,batch):
    target=batch
    chromosome_file=open(f'{outputFile}.chromosomes','w')
    chromosome_file.close()
    counter=0
    tmp_reads=open(f'{outputFile}.tmp','w')
    with open(samFile,'r') as f:
        for line in f:
            if not line.startswith('@'):
                counter+=1
                if counter%target==0:
                    print('\t\tprocessing read', counter,' '*20,end='\r')
                    tmp_reads.close()
                    os.system(f'python3 {MandoPath}/utils/emtreyProcessSamBatch.py -t {threads} -o {outputFile} {mando}')
                    tmp_reads=open(f'{outputFile}.tmp','w')
                tmp_reads.write(line)
            else:
                if line.startswith('@SQ'):
                    a=line.rstrip().split('\t')
                    chr=a[1].split(':')[1]
                    length=int(a[2].split(':')[1])
                    chromosome_file=open(f'{outputFile}.chromosomes','a')
                    chromosome_file.write(f'{chr}\t{length}\n')
                    chromosome_file.close()
    tmp_reads.close()
    os.system(f'python3 {MandoPath}/utils/emtreyProcessSamBatch.py -t {threads} -o {outputFile} {mando}')
    os.system(f'rm {outputFile}.tmp')


def main():
    delegatingPslConversion(inFile,batch)

main()
