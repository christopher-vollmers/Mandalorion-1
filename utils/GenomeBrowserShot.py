import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import matplotlib.patches as mplpatches
from matplotlib.collections import PatchCollection

viridis5 = (253/255, 231/255, 37/255)
viridis4 = (94/255, 201/255, 98/255)
viridis3 = (33/255, 145/255, 140/255)
viridis2 = (59/255, 82/255, 139/255)
viridis1 = (68/255, 1/255, 84/255)

R1=np.linspace(viridis1[0],viridis2[0],26)[:-1]
G1=np.linspace(viridis1[1],viridis2[1],26)[:-1]
B1=np.linspace(viridis1[2],viridis2[2],26)[:-1]

R2=np.linspace(viridis2[0],viridis3[0],26)[:-1]
G2=np.linspace(viridis2[1],viridis3[1],26)[:-1]
B2=np.linspace(viridis2[2],viridis3[2],26)[:-1]

R3=np.linspace(viridis3[0],viridis4[0],26)[:-1]
G3=np.linspace(viridis3[1],viridis4[1],26)[:-1]
B3=np.linspace(viridis3[2],viridis4[2],26)[:-1]

R4=np.linspace(viridis4[0],viridis5[0],26)
G4=np.linspace(viridis4[1],viridis5[1],26)
B4=np.linspace(viridis4[2],viridis5[2],26)

R=np.concatenate((R1,R2,R3,R4),axis=None)
G=np.concatenate((G1,G2,G3,G4),axis=None)
B=np.concatenate((B1,B2,B3,B4),axis=None)

RB1=(225/255,13/255,50/255)
RB2=(242/255,50/255,54/255)
RB3=(239/255,99/255,59/255)
RB4=(244/255,138/255,30/255)
RB5=(248/255,177/255,61/255)
RB6=(143/255,138/255,86/255)
RB7=(32/255,100/255,113/255)
RB8=(42/255,88/255,132/255)
RB9=(56/255,66/255,156/255)
RB10=(84/255,60/255,135/255)
RB11=(110/255,57/255,115/255)
RB12=(155/255,42/255,90/255)


def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='plots stuff',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--transcriptModels','-t', type=str, action='store',
                          help='gtf files for bottom panel')
    parser.add_argument('--transcriptRpm','-f', type=str, action='store',
                          help='Mando .rpm file for bottom panel')
    parser.add_argument('--genomeAnnotation', '-g', type=str, action='store',
                          help='gtf file for second panel')
    parser.add_argument('--normalization', '-n', type=str, action='store', default='all',choices=['isoform','gene','geneSample','all'],
                          help='set how heatmap should be normalized')
    parser.add_argument('--introns', '-i', action='store_true',
                          help='set if you want to include introns')
    parser.add_argument('--range', '-r', type=str, action='store',
                        help='''Comma separated coordinates for genome window to show: chromosome,start,end''')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()



plt.style.use('BME163')




def read_psl(psl,genomeRange,color1,color2,direction_selection):
    targetChrom,targetStart,targetEnd=genomeRange
    colors={}
    if not color1:
        colors['+']=RB3
        colors['-']=RB6
    else:
        colors['+']=color1
        colors['-']=color2

    reads=[]
    with open(psl) as f:
        for line in f:
            use=False
            splitLine=line.split('\t')
            chromosome,start,end,direction,name=splitLine[13],int(splitLine[15]),int(splitLine[16]),splitLine[8],splitLine[9]
            if chromosome==targetChrom:
                if targetStart<start<targetEnd or targetStart<end<targetEnd:
                    use=True
                elif targetStart>start and end>targetEnd:
                    use=True
            if use:
                blockSizes=np.array(splitLine[18].split(',')[:-1],dtype=int)
                blockStarts=np.array(splitLine[20].split(',')[:-1],dtype=int)
                blockHeights=[0.75]*len(blockStarts)
                color=colors[direction]
                use1=True
                if direction_selection:
                    use1=False
                    if direction == direction_selection:
                        use1=True
                if use1:
                    reads.append([chromosome,start,end,name,direction,blockStarts,blockSizes,blockHeights,color])
    return reads

def read_gtf(gtf,genomeRange,color1,color2,direction_selection):
    targetChrom,targetStart,targetEnd=genomeRange
    colors={}
    if not color1:
        colors['+']=RB1
        colors['-']=RB6
    else:
        colors['+']=color1
        colors['-']=color2

    heightDict={}
    heightDict['exon']=0.375
    heightDict['CDS']=0.75

    transcript_dict={}
    with open(gtf) as f:
        for line in f:
            splitLine=line.strip().split('\t')
            if len(splitLine)>7:
                if splitLine[0]==targetChrom:
                    if splitLine[2] in ['exon','CDS']:
                        transcriptID=splitLine[8].split('transcript_id "')[1].split('"')[0]
                        if transcriptID not in transcript_dict:
                            transcript_dict[transcriptID]=[]
                        transcript_dict[transcriptID].append((splitLine[2],int(splitLine[3]),int(splitLine[4]),splitLine[6]))
    reads=[]
    for transcript_id,entries in transcript_dict.items():
        starts=[]
        ends=[]
        for entry in entries:
            starts.append(entry[1])
            ends.append(entry[2])
        start=min(starts)
        end=max(ends)
        use=False
        if targetStart<start<targetEnd or targetStart<end<targetEnd:
            use=True
        elif targetStart>start and end>targetEnd:
            use=True
        if use:
            blockStarts=[]
            blockSizes=[]
            blockHeights=[]
            for entry in transcript_dict[transcript_id]:
                type=entry[0]
                direction=entry[3]
                blockStarts.append(int(entry[1]))
                blockSizes.append(int(entry[2])-int(entry[1]))
                if heightDict:
                    blockHeight=heightDict[type]
                else:
                    blockHeight=0.5
                blockHeights.append(blockHeight)
                color=colors[direction]
            use1=True
            if direction_selection:
                use1=False
                if direction == direction_selection:
                    use1=True
            if use1:
                reads.append([targetChrom,start,end,transcript_id,direction,blockStarts,blockSizes,blockHeights,color])

    return reads

def nostack(reads,bottom):
    level=bottom
    for read in reads:
        top=level
        read.append(level)
        level+=1
    return reads,top+1

def stack(reads,bottom):
    top=bottom
    for level in range(bottom,bottom+len(reads),1):
        previous=-1
        for read in reads:
            if len(read)==9:

                top=level
                start=read[1]
                end=read[2]

                if start>previous:
                    read.append(level)
                    previous=end

    return reads,top+1

def plot_reads(reads,reverse,bottom,subsample,coverage,left,right,stacking):
    bottom+=1
    if subsample:
        subReads=[]
        indeces = np.random.choice(np.arange(0, len(reads)),
                                   min(len(reads), subsample), replace=False)
        for index in indeces:
            read = reads[index]
            subReads.append(read)
    reads=subReads

    if not stacking:
        if reverse:
            reads=sorted(reads,key=lambda x:(x[4],x[2]),reverse=False)
        else:
            reads=sorted(reads,key=lambda x:(x[4],x[1]),reverse=False)
    else:
        if reverse:
            reads=sorted(reads,key=lambda x:(x[2]),reverse=False)
        else:
            reads=sorted(reads,key=lambda x:(x[1]),reverse=False)


    rectangles=[]
    if stacking:
        stackedReads,top=stack(reads,bottom)
    else:
        stackedReads,top=nostack(reads,bottom)
    for read in stackedReads:
        targetChrom,start,end,type,direction,blockStarts,blockSizes,blockHeights,color,y_pos=read
        rectangles.append(mplpatches.Rectangle((start,y_pos-0.025),end-start,0.05,facecolor=color,edgecolor='black',linewidth=0.1))
        for pos in range(len(blockStarts)):
            rectangles.append(mplpatches.Rectangle((blockStarts[pos],y_pos-(blockHeights[pos]/2)),blockSizes[pos],blockHeights[pos],facecolor=color,edgecolor='black',linewidth=0.1))
            for pos2 in range(blockStarts[pos],blockStarts[pos]+blockSizes[pos],1):
                coverage.add(pos2)
#    print(stackedReads,top)

    rectangles.append(mplpatches.Rectangle((left,top),right-left,0.1,facecolor='black',edgecolor='black',linewidth=0.1))

    return rectangles,top,coverage,stackedReads



# left,bottom, width,height

def compile(coverage,genomeRange):

    on=0
    off=0
    areas=[]
    targetChrom,targetStart,targetEnd=genomeRange
    areaStart,areaEnd=False,False
    for pos in np.arange(targetStart,targetEnd,1):
        if pos in coverage:

            if on:
                areaEnd=pos+45
            else:

                areaStart=pos-45
                on=100
        else:
            if on:
                on-=1
            if not on:
                if areaStart and areaEnd:
                    areas.append((areaStart,areaEnd))
                    areaStart,areaEnd=False,False

    if on:
        areas.append((areaStart,areaEnd))
    return areas

def make_panels(areas,rpH,panelBottom):
    overall=0
    panels={}
    for start,end in areas:
        overall+=(end-start)
    previous=0.02
    for start,end in areas:
        fraction=((end-start)/overall)*0.70
        panels[(start,end)]=[previous,panelBottom,fraction,rpH]
        previous+=(fraction+min(0.005,0.1/len(areas)))
        kb = (1000/overall)*0.7
    return panels, kb

def draw_borders(panel1,start,end,limits,top1):
    left=min(limits)
    right=max(limits)
    bottom=0.5
    top1-=0.5

    if left!=start:
        panel1.spines['left'].set_visible(False)
        panel1.plot([start,start],[bottom,top1],color='black',linestyle='--',lw=0.25)
    if right!=end:
        panel1.spines['right'].set_visible(False)
        panel1.plot([end,end],[bottom,top1],color='black',linestyle='--',lw=0.25)

#    panel1.plot([left,left],[bottom,top1],color='black',lw=1)
#    panel1.plot([right,right],[bottom,top1],color='black',lw=1)
#    panel1.plot([left,right],[top1,top1],color='black',lw=1)
#    panel1.plot([left,right],[bottom,bottom],color='black',lw=1)
    panel1.set_xlim(start,end)
    panel1.set_ylim(bottom,top1)

def plot_heatmap(panel,reads,frac1,norm):
    rectangles=[]
    yDict={}
    yList=[]
    for read in reads:
        name, y = read[3],read[9]
        yDict[name]=y
        yList.append(y)
    geneFracs={}
    geneList={}
    with open(frac1) as f:
        info=f.readline().strip()
        samples=info.split('\t')[2:]
        totalFracs=[]
        for line in f:
            splitLine = line.strip().split('\t')
            name,gene,fracs = splitLine[0],splitLine[1],np.array(splitLine[2:],dtype=float)
            if name in yDict:
                if gene not in geneFracs:
                    geneFracs[gene]=np.zeros(len(samples))
                    geneList[gene]=[]
                totalFracs+=list(fracs)
                geneFracs[gene]+=fracs
                geneList[gene]+=list(fracs)
    with open(frac1) as f:
        info=f.readline().strip()
        samples=info.split('\t')[2:]
        for line in f:
            splitLine = line.strip().split('\t')
            name,gene,fracs = splitLine[0],splitLine[1],np.array(splitLine[2:],dtype=float)
            if name in yDict:
                y_pos=yDict[name]
                x_pos=-1
                for index in range(len(fracs)):
                    frac=fracs[index]
                    x_pos+=1
                    if norm == 'isoform':
                        maximum = np.max(fracs)
                    if norm == 'all':
                        maximum = np.max(totalFracs)
                    if norm == 'gene':
                        maximum = np.max(geneList[gene])
                    if norm == 'geneSample':
                        maximum = geneFracs[gene][index]
#                        step=int(((float(frac)-minimum)/(maximum-minimum))*100)
                    if maximum>0:
                        step=int((float(frac)/maximum)*100)
                        color=R[step],G[step],B[step]
                    else:
                        color='grey'
                    rectangles.append(mplpatches.Rectangle((x_pos,y_pos-0.5),1,1,facecolor=color,edgecolor='black',linewidth=0.3))
                    #panel.text(x_pos+0.5,y_pos,str(count),va='center',ha='center',fontsize=4)
    patches = PatchCollection(rectangles,match_original=True)
    panel.add_collection(patches)
    trunc_samples=[]

    x=0
    for sample in samples:
        x+=1
        trunc_samples.append(sample.split('/')[-1])
#        panel.plot([x,x],[0,max(yList)+1],linewidth=0.5,color='white')
    panel.set_xlim(0,len(samples))
    panel.set_ylim(0.5,max(yList)+0.5)
    panel.set_xticks(np.arange(0.5,len(samples)+0.5,1),trunc_samples,rotation=90,fontsize=4)
    panel.tick_params(bottom=False, right=False,
                           left=False,top=True,
                           labelbottom=False,labelleft=False,
                           labeltop=True,labelright=False)
def main():
    parser=parse_args()
    figureWidth=4
    figureHeight=4
    panelHeight=1
    panelHeight2=0.5
    plt.figure(figsize=(figureWidth,figureHeight))


    group1,group2,genomeRange,frac1,normalization,introns=parser.transcriptModels,parser.genomeAnnotation,parser.range,parser.transcriptRpm,parser.normalization,parser.introns
    groups1=False
    if ',' in group1:
        groups1=group1.split(',')
    else:
        groups1=[group1]

    groups2=False
    if ',' in group2:
        groups2=group2.split(',')
    else:
        groups2=[group2]


    values=genomeRange.split(',')
    genomeRange=(values[0],int(values[1]),int(values[2]))
    print(genomeRange)

    group1Reads=[]
    for group1 in groups1:
        if group1.endswith('.gtf'):
            group1Reads.append(read_gtf(group1,genomeRange,RB10,RB5,False))
        elif group1.endswith('.psl'):
            group1Reads.append(read_psl(group1,genomeRange,RB10,RB5,False))

    group2Reads=[]
    for group2 in groups2:
        if group2.endswith('.gtf'):
            group2Reads.append(read_gtf(group2,genomeRange,RB10,RB5,False))
        elif group2.endswith('.psl'):
            group2Reads.append(read_psl(group2,genomeRange,RB10,RB5,False))


    coverage=set()
    top1=0
    rectangles=[]
    for reads in group1Reads:
        rectangle,top1,coverage,reads1=plot_reads(reads,False,top1,200,coverage,genomeRange[1],genomeRange[2],False)
        rectangles.append(rectangle)

    top2=0
    rectangles2=[]
    for reads in group2Reads:
        rectangle,top2,coverage,reads2=plot_reads(reads,False,top2,200,coverage,genomeRange[1],genomeRange[2],True)
        rectangles2.append(rectangle)


    if introns:
        for i in range(genomeRange[1],genomeRange[2],1):
            coverage.add(i)

    areas=compile(coverage,genomeRange)

    panels,kb=make_panels(areas,panelHeight/figureHeight,0.1)
    panels2,kb=make_panels(areas,panelHeight2/figureHeight,0.36)

    limits=[]
    for (start,end),panel_info in panels.items():
        limits.append(start)
        limits.append(end)

    right_most = []
    for (start,end),panel_info in panels2.items():
        panel2=plt.axes(panel_info,frameon=True)
        right_most.append(panel_info[0]+panel_info[2])

        for index in np.arange(0,len(groups2),1):
            rectangle = rectangles2[index]
            name=groups2[index].split('/')[-1]
            patches = PatchCollection(rectangle,match_original=True)
            panel2.add_collection(patches)
        draw_borders(panel2,start,end,limits,top2)
        panel2.tick_params(bottom=False, right=False,
                           left=False,top=False,
                           labelbottom=False,labelleft=False,
                           labeltop=False,labelright=False)

    for (start,end),panel_info in panels.items():
        panel1=plt.axes(panel_info,frameon=True)

        for index in np.arange(0,len(groups1),1):
            rectangle = rectangles[index]
            name=groups1[index].split('/')[-1]
            patches = PatchCollection(rectangle,match_original=True)
            panel1.add_collection(patches)
        draw_borders(panel1,start,end,limits,top1)
        panel1.tick_params(bottom=True, right=False,
                           left=False,top=False,
                           labelbottom=True,labelleft=False,
                           labeltop=False,labelright=False)
#        panel1.set_xticks([start,end],[f'\n{start:,}',f'{end:,}\n'],rotation=90,fontsize=4)
        panel1.set_xticks([start],[f'{start:,}'],rotation=90,fontsize=4)
    panel3=plt.axes([max(right_most)+0.01,0.1,0.98-(max(right_most)+0.02),panelHeight/figureHeight],frameon=True)
    plot_heatmap(panel3,reads1,frac1,normalization)
    panel4=plt.axes([0.02,0.47,kb,0.05],frameon=False)
    panel4.plot([0,1],[1,1],linewidth=1,color='black')
    panel4.text(0.5,1,'1kb',va='bottom',ha='center')
    panel4.set_xlim(0,1)
    panel4.set_ylim(0,2)
    panel4.tick_params(bottom=False, right=False,
                       left=False,top=False,
                       labelbottom=False,labelleft=False,
                       labeltop=False,labelright=False)

    plt.savefig('Final_assignment.png',dpi=2400)
main()
