#! /usr/bin/env python3

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse
from decimal import Decimal
import numpy
import random
from signal import signal, SIGPIPE, SIG_DFL
import os.path
import json


signal(SIGPIPE, SIG_DFL) # Handle broken pipes


version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

def main():

    parser = argparse.ArgumentParser(description="Create random SNVs, indels and CNVs for each subclone in a tumour sample.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(version['__version__']))
    parser.add_argument('-j', '--json', dest='jsonfile', required=True, type=str, help='Json file with parameters')

    args = parser.parse_args()
    with open(args.jsonfile,'r') as file:
        parameters=json.load(file)

    #set default parameter values and give error/warning messages
    if "prefix" not in parameters:
        parameters['prefix'] = ''
    if "reference" not in parameters:
        print('Error: No input genome fasta file provided.')
    if os.path.exists(parameters['reference'] + '.fai'):
        parameters['fai']=(parameters['reference'] + '.fai')
    else:
        print('Error: No fai index file for genome.')       
    if "directory" not in parameters:
        print('Warning: No output directory given, using current directory.')
        parameters['directory']='./'
    if "structure" not in parameters:
        print('Warning: No tumour structure given, using "clone1,0.2,germline,clone2,0.8,clone1"')
        parameters['structure']="clone1,0.2,germline,clone2,0.8,clone1"
    if "snvgermline" not in parameters:
        parameters['snvgermline']=0.00014
        print('Warning: No germline SNV rate given, using '+str(parameters['snvgermline'])+'.')
    if "indgermline" not in parameters:
        parameters['indgermline']=0.000014
        print('Warning: No germline InDel rate given, using '+str(parameters['indgermline'])+'.')
    if "cnvrepgermline" not in parameters:
        parameters['cnvrepgermline']=160
        print('Warning: No germline CNV number given, using '+str(parameters['cnvrepgermline'])+'.')
    if "cnvdelgermline" not in parameters:
        parameters['cnvdelgermline']=1000
        print('Warning: No germline CNV number given, using '+str(parameters['cnvdelgermline'])+'.')
    if "aneuploid" not in parameters:
        parameters['aneuploid']=2
        print('Warning: No somatic aneuploid events number given, using '+str(parameters['aneuploid'])+'.')
    if "snvsomatic" not in parameters:
        parameters['snvsomatic']=0.00001
        print('Warning: No somatic SNV rate given, using '+str(parameters['snvsomatic'])+'.')
    if "indsomatic" not in parameters:
        parameters['indsomatic']=0.000001
        print('Warning: No somatic InDel rate given, using '+str(parameters['indsomatic'])+'.')
    if "cnvrepsomatic" not in parameters:
        parameters['cnvrepsomatic']=250
        print('Warning: No somatic replication CNV number given, using '+str(parameters['cnvrepsomatic'])+'.')
    if "cnvdelsomatic" not in parameters:
        parameters['cnvdelsomatic']=250
        print('Warning: No somatic deletion CNV number given, using '+str(parameters['cnvdelsomatic'])+'.')
    if "cnvgermlinemean" not in parameters:
        parameters['cnvgermlinemean']=-10
        print('Warning: No lognormal mean for germline CNV lengths given, using '+str(parameters['cnvgermlinemean'])+'.')
    if "cnvgermlinevariance" not in parameters:
        parameters['cnvgermlinevariance']=3
        print('Warning: No lognormal variance for germline CNV lengths given, using '+str(parameters['cnvgermlinevariance'])+'.')
    if "cnvgermlinemultiply" not in parameters:
        parameters['cnvgermlinemultiply']=1000000
        print('Warning: No multipication factor for germline CNV length given, using '+str(parameters['cnvgermlinemultiply'])+'.')
    if "cnvsomaticmean" not in parameters:
        parameters['cnvsomaticmean']=-1
        print('Warning: No lognormal mean for somatic CNV lengths given, using '+str(parameters['cnvsomaticmean'])+'.')
    if "cnvsomaticvariance" not in parameters:
        parameters['cnvsomaticvariance']=3
        print('Warning: No lognormal variance for somatic CNV lengths given, using '+str(parameters['cnvsomaticvariance'])+'.')
    if "cnvsomaticmultiply" not in parameters:
        parameters['cnvsomaticmultiply']=1000000
        print('Warning: No multipication factor for somatic CNV length given, using '+str(parameters['cnvsomaticmultiply'])+'.')
    if "indmean" not in parameters:
        parameters['indmean']=-2
        print('Warning: No mean for InDel lengths given, using '+str(parameters['indmean'])+'.')
    if "indvariance" not in parameters:
        parameters['indvariance']=2
        print('Warning: No variance for InDel lengths given, using '+str(parameters['indvariance'])+'.')
    if "indmultiply" not in parameters:
        parameters['indmultiply']=1
        print('Warning: No multipication factor for InDel length given, using '+str(parameters['indmultiply'])+'.')
    if "cnvcopiesmean" not in parameters:
        parameters['cnvcopiesmean']=1
        print('Warning: No mean cnv copy number given, using '+str(parameters['cnvcopiesmean'])+'.')
    if "cnvcopiesvariance" not in parameters:
        parameters['cnvcopiesvariance']=0.5
        print('Warning: No variance for cnv copy number given, using '+str(parameters['cnvcopiesvariance'])+'.')
    if "chromosomes" not in parameters:
        parameters['chromosomes']='all'
        print('Warning: No chromosomes given, using '+str(parameters['chromosomes'])+'.')
    if "dbsnp" not in parameters:
        print('Warning: No dnSNP file given, all variants will be randomly generated.')
        parameters['dbsnp']='none'
        parameters['dbsnpindelproportion']=0
        parameters['dbsnpsnvproportion']=0
    else:
        if "dbsnpsnvproportion" not in parameters:
            parameters['dbsnpsnvproportion']=0.9
            print('Warning: No proportion given for germline SNVs taken from dbSNP, using '+str(parameters['dbsnpsnvproportion'])+'.')
        if "dbsnpindelproportion" not in parameters:
            parameters['dbsnpindelproportion']=0.5
            print('Warning: No proportion given for germline InDels taken from dbSNP, using '+str(parameters['dbsnpindelproportion'])+'.')


    #Functions for reading in data----------------------------------------------------------------------------------

    def readinclones(parameters):    #formats structure parameter into dictionary
        if os.path.exists(parameters['structure']):
            with open(parameters['structure'],'r') as file:
                structure=file.readline().strip().split(',')
        else:
            structure=parameters['structure'].split(',')

        clones={}
        for i in range(1,int((len(structure)/3)+1)):
            clones[structure[i*3-3]]=[structure[i*3-2],structure[i*3-1]]
        return(clones)

    def readindbsnp(dbsnp,gen):    #formats structure parameter into dictionary
        with open(dbsnp,'r') as file:
            dbsnvs={}
            dbindels={}
            for chro in gen:
                dbsnvs[chro]=[]
                dbindels[chro]=[]
            for line in file:
                if line.startswith('#'):
                    continue
                l=line.strip().split("\t")
                if l[0] in gen:
                    ref = l[3].split(',')
                    alt = l[4].split(',')
                    for i in ref:
                        for ii in alt:
                            if len(ii)==1 and len(i)==1:    #if its an snv
                                dbsnvs[l[0]].append([l[1],i,ii])
                            elif len(ii)==1 or len(i)==1:   #if its an indel, but not multiple bases being replaced by multiple bases (too much hassel)
                                dbindels[l[0]].append([l[1],i,ii])
        return(dbsnvs,dbindels)

    def readinfai(chromosomes,fai,referencefile): #reads in reference and fai files for required chromosomes into dictionaries
        if chromosomes==['all']:
            keepchromos=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
        else:
            keepchromos=chromosomes
        with open(fai,'r') as file:
                gen = dict([(line.strip().split("\t")[0], float(line.strip().split("\t")[1])) for line in file])
        genkeys=list(gen.keys())
        for chro in genkeys:
            if chro not in (keepchromos):
                del gen[chro]
        reference={}
        reflist=[]
        chromo=''
        with open(referencefile,'r') as ref:
            for l in ref:
                if l.startswith('>'):
                    if chromo!='': reference[chromo]=''.join(reflist)
                    reflist=[]
                    chromo=l.strip().split(" ")[0][1:]
                    reference[l.strip().split(" ")[0]]=''
                else:
                    reflist.append(l.strip())
            reference[chromo]=''.join(reflist)
        reflist=''
        for chro in list(reference.keys()):
            if chro not in (keepchromos): del reference[chro]
        return(gen,reference)

    #Functions for variant generation-------------------------------------------------------------------------------------------------

    def choosechromosome(gen,chrohaps):
        totlength=sum(gen[i]*len(chrohaps[i]) for i in gen.keys()) #get total length of all chromosome copies
        if totlength==0:
            print('All chromosomes have been deleted, no genome left!')
            exit()
        probs=[float(gen[i]*len(chrohaps[i]))/totlength for i in list(gen.keys())] #get probability of each type of chromosome being choosen
        chro=numpy.random.choice(list(gen.keys()), p=probs)
        hap=numpy.random.choice(chrohaps[chro])
        return(chro,hap)

    def createaneu (gen,chrohaps):  #create information for an aneuploid variant
        chro=numpy.random.choice(list(gen.keys())) #choose chromosome randomly with even probabilities
        hap=numpy.random.choice(chrohaps[chro])
        copy=int(numpy.random.choice([0,2,3]))   #choose copy number of chromosome
        print(['aneu',chro,hap,copy])
        return ['aneu',chro,hap,copy]

    def createcnv(gen,chrohaps,vartype,somorger):   #create information for a cnv variant
        chro,hap=choosechromosome(gen,chrohaps)
        startbase='N'
        endbase='N'
        while (startbase=='N' or startbase=='n') and (endbase=='N' or endbase=='n'):
            length=0
            if somorger=='germline':
                while length<51 or length>gen[chro]/2: length=int(Decimal(numpy.random.lognormal(parameters['cnvgermlinemean'],parameters['cnvgermlinevariance'],1)[0])*parameters['cnvgermlinemultiply'])
            else:
                while length<51 or length>gen[chro]/2: length=int(Decimal(numpy.random.lognormal(parameters['cnvsomaticmean'],parameters['cnvsomaticvariance'],1)[0])*parameters['cnvsomaticmultiply'])
            position=random.randint(1,gen[chro]-length+1)
            startbase=reference[chro][position-1]
            endbase=reference[chro][position-1+length-1]
        copy=0
        if vartype=='cnvrep':
            while copy==0 or copy==1: copy=int(numpy.random.lognormal(parameters['cnvcopiesmean'],parameters['cnvcopiesvariance'],1))
        return ['cnv',chro,hap,position,length,copy]

    def createsnv(gen,chrohaps,dbsnvs,pro):    #create information for an snv variant
        chro,hap=choosechromosome(gen,chrohaps)
        if pro=='':
            i=False
        else:
            i=numpy.random.choice([True,False],1,p=[float(pro),(1-float(pro))])[0]
        if i == True:   #if taking variant form dnindels
            l=dbsnvs[chro][numpy.random.choice(range(len(dbsnvs[chro])))]
            position=int(l[0])
            ref=l[1]
            alt=l[2]
        else:
            ref='N'
            while ref=='n' or ref=='N':
                position=random.randint(1,gen[chro])
                ref=reference[chro][position-1]
            substitutions={}
            substitutions['A']=['T','G','C']
            substitutions['T']=['A','G','C']
            substitutions['C']=['T','G','A']
            substitutions['G']=['T','A','C']
            substitutions['a']=['t','g','c']
            substitutions['t']=['a','g','c']
            substitutions['c']=['t','g','a']
            substitutions['g']=['t','a','c']
            alt=str(numpy.random.choice(substitutions[ref]))
        return ['snv',chro,hap,position,ref,alt]

    def createind(gen,chrohaps,dbindels,pro):    #create information for an indel variant
        chro,hap=choosechromosome(gen,chrohaps)
        if pro=='':
            i=False
        else:
            i=numpy.random.choice([True,False],1,p=[float(pro),(1-float(pro))])[0]
        if i == True:   #if taking variant form dnindels
            l=dbindels[chro][numpy.random.choice(range(len(dbindels[chro])))]
            print(l)
            position=int(l[0])
            ref=l[1]
            alt=l[2]
            if len(ref) != 1:
                length=len(ref)
                iod='d'
            else:
                length=len(alt)
                iod='i'
        else:
            length=101
            while length>50:
                length=int(round((float(numpy.random.lognormal(parameters['indmean'],parameters['indvariance'],1)[0])*parameters['indmultiply'])+0.5,0))
            seq=50*'N'
            refplusone='N'
            while (seq.count('N')+seq.count('n')>float(length)/4) or refplusone=='n' or refplusone=='N':    #keep getting an indel until it doesn't contain more than 1/4 'N's or base after position isn't an N (indels start on base after position)
                position=random.randint(1,gen[chro]-length)
                refplusone=reference[chro][position-1+1]
                iod=numpy.random.choice(['i','d'])
                if iod=='i':
                    ref=reference[chro][position-1]
                    seqchro,dontneed=choosechromosome(gen,chrohaps)
                    seqpos=random.randint(1,gen[seqchro]-length+1)
                    seq=ref + reference[seqchro][seqpos-1:seqpos+length-1]  #get sequence from elsewhere in the genome
                    alt=seq
                    ref=alt[0]
                else:
                    seq=reference[chro][position-1:position+length]#get deleted sequence
                    ref=seq
                    alt=ref[0]
        print(['indel',chro,hap,position,length,ref,alt,iod])
        return ['indel',chro,hap,position,length,ref,alt,iod]

    def getcnv(gen,lists,vartype,somorger):
        keep=False
        c=1
        while keep==False:  #keep getting variant until it fits
            v=createcnv(gen,lists[1],vartype,somorger)
            keep=True
            for x in lists[3][v[1]+v[2]]:    #for each breakpoint pair in cnv dictionary
                if v[3] < x[0] and  v[3]+v[4]-1 >= x[0] and v[3]+v[4]-1 < x[1]:   #if start position out but end position in. (-1 is added to length as the position base is included in the length)
                    keep=False
                    c+=1
                    break
                elif v[3]+v[4]-1 > x[1] and  v[3] <= x[1] and v[3] > x[0]:   #if end position out but start position in
                    keep=False
                    c+=1
                    break
            for x in lists[4][v[1]+v[2]]:    #for each breakpoint pair in deletions dictionary
                if (v[3] >= x[0] and  v[3] <= x[1]) or (v[3]+v[4]-1 >= x[0] and v[3]+v[4]-1 <= x[1]):   #if start or end positions in deleted region
                    keep=False
                    c+=1
                    break
                elif v[3]==x[0]: #if cnv starts at same base as existing cnv
                    keep=False
                    c+=1
                    break
            if c>100:
                print('Warning, not enough room in genome for so many CNVs. Program may not end. Kill me now...')

        lists[0].append(v)   #add variant to variants list
        if v[5]!=0: #if cnv is not a deletion
            lists[3][v[1]+v[2]].append([v[3],v[3]+v[4]-1]) #add to the cnv breakpoints dictionary for key(chromosome+haplotype), [start,end]
        if v[5]==0: #if cnv is a deletion
            lists[4][v[1]+v[2]].append([v[3],v[3]+v[4]-1]) #add to the deleted regions dictionary for key(chromosome+haplotype), [start,end]
        return(lists)

    def getind(gen,lists,dbindels,dbsnpindelproportion):
        #when refering to the location of an indel, +1 is used as the position refers to the previous base
        #-1 is added to length as the first base of the indel (position+1) is included in the length
        keep=False
        while keep==False:
            v=createind(gen,lists[1],dbindels,dbsnpindelproportion)
            keep=True
            if v[3] in lists[2][v[1]+v[2]]:  #if position exists in snv/indel dictionary
                    keep=False
            if v[7]=='i':
                for x in lists[4][v[1]+v[2]]:    #for each breakpoint pair in deletions dictionary
                    if v[3] <= x[1] and  v[3]+1 >= x[0]:   #if start position (or previous base) in deleted region
                        keep=False
            if v[7]=='d':
                for x in lists[4][v[1]+v[2]]:    #for each breakpoint pair in deletions dictionary
                    if (v[3] <= x[1] and  v[3]+1 >= x[0]) or (v[3]+1 +v[4]-1 <= x[1] and  v[3]+1+v[4]-1 >= x[0]):   #if start position(or previous base) or end position in deleted region
                        keep=False
                for x in lists[3][v[1]+v[2]]:    #for each breakpoint pair in cnv breakpoints dictionary
                    if (v[3] <= x[0] and  v[3]+1 +v[4]-1  >= x[0]) or (v[3] <= x[1] and  v[3]+1 +v[4]-1 >= x[1]):   #if cnv start position or end position in deleted region
                        keep=False
                for x in range(v[3]+1,v[3]+1 +v[4]-1):  #if indel covers an existing snv
                    if x in lists[2][v[1]+v[2]]:
                        keep=False

        lists[0].append(v)   #add variant to variants list
        lists[2][v[1]+v[2]].append(v[3]) #add to the indel dictionary for key(chromosome+haplotype)
        if v[6]=='d': #if indel is a deletion
            lists[4][v[1]+v[2]].append([v[3]+1,v[3]+1+v[4]-1]) #add to the deleted regions dictionary for key(chromosome+haplotype), [start,end] - previous base is also added to simplify the genome writing stage (don't want a cnv starting/ending between previous base and deleted region)
        return(lists)

    def getsnv(gen,lists,dbsnvs,dbsnpsnvproportion):
        keep=False
        while keep==False:
            v=createsnv(gen,lists[1],dbsnvs,dbsnpsnvproportion)
            keep=True
            if v[3] in lists[2][v[1]+v[2]]:  #if position exists in dictionary
                keep=False
            for x in lists[4][v[1]+v[2]]:    #for each breakpoint pair in dictionary
                if v[3] <= x[1] and  v[3] >= x[0]:   #if position in deleted region
                    keep=False
        lists[0].append(v)   #add variant to variants list
        lists[2][v[1]+v[2]].append(v[3])
        return(lists)

    def getaneu(gen,lists):
        v=createaneu(gen,lists[1])
        lists[0].append(v)
        #update chromosome haplotypes
        newhaps=[]
        for i in range(1,v[3]+1):
            newhaps.append(v[2]+'-'+str(i))
        for hap in newhaps:
            lists[1][v[1]].append(hap)   #add haplotype for each copy number - eg. A,B becomes A,B,A1,A2,A3
        lists[1][v[1]].remove(v[2])       #remove original
        #update lists to replicate for each new haplotype
        for l in [2,3,4]:
            for chrohap in list(lists[l].keys()):
                if chrohap==v[1]+v[2]:
                    for hap in newhaps:
                        lists[l][v[1]+hap]=lists[l][v[1]+v[2]]
        #add empty keys for each haplotype
            for hap in lists[1][v[1]]:
                if v[1]+hap not in lists[l]:
                    lists[l][v[1]+hap]={}
        if v[3]==0:
            lists[4][v[1]+v[2]]=[1,gen[v[1]]]
        return(lists)

    def chooseclone(clones,totdis):
        clo=numpy.random.choice(list(clones.keys()),p=[float(clones[c][0])/totdis for c in clones])
        return(clo)

    def shufflelist(shuflist):
        shuflist=numpy.random.choice(shuflist,replace=False,size=len(shuflist))
        return(shuflist)


    #Functions for writing output files--------------------------------------------------------------------------------

    def writevariantfile(directory,prefix,lists,clo):
        with open(directory + '/' + prefix + clo+'variants.txt','w+') as file:
            i=0
            for v in lists[0]:
                i+=1
                for x in v:
                    file.write(str(x) + '\t')
                file.write(str(i) + '\n')

    def writeaneuploidfile(directory,prefix,lists,clo,germorsom,gen):
        with open(directory + '/' + prefix + clo + 'aneuploid.txt','w+') as file:
            for v in lists[0]:
                if v[0]=='aneu':
                    file.write(v[1]+'\t'+v[2]+'\t'+str(v[3])+'\n')
            file.write('Total chromosomes:' + '\n')
            for chromo in gen:
                for hap in lists[1][chromo]:
                    file.write(chromo+hap + '\n')


    #Preparation---------------------------------------------------------------------------------------------

    #Read in clones and reference genome
    clones=readinclones(parameters)    #get dictionary of specified clones
    if type(parameters['chromosomes'])==list:
        chromosomes=parameters['chromosomes']
    else:
        chromosomes=[parameters['chromosomes']]
    gen,reference=readinfai(chromosomes,parameters['fai'],parameters['reference'])  #get dictionaries of genome lengths and sequences
    print('Tumour clones : ',clones)

    #read in dbsnp if given
    if parameters['dbsnp'] != 'none':
        dbsnvs,dbindels=readindbsnp(parameters['dbsnp'],gen)
    else:
        dbsnvs={}
        dbindels={}

    #Get total number of each variant type for somatic and germline genomes
    snvgernum=round(sum(list(gen.values()))*parameters['snvgermline'])
    snvsomnum=round(sum(list(gen.values()))*parameters['snvsomatic'])
    indgernum=round(sum(list(gen.values()))*parameters['indgermline'])
    indsomnum=round(sum(list(gen.values()))*parameters['indgermline'])

    print('Number of germline SNVs : ',snvgernum)
    print('Number of somatic SNVs : ',snvsomnum)
    print('Number of germline indels : ',indgernum)
    print('Number of somatic indels : ',indsomnum)
    print('Number of germline replication CNVs : ',parameters['cnvrepgermline'])
    print('Number of somatic replication CNVs : ',parameters['cnvrepsomatic'])
    print('Number of germline deletion CNVs : ',parameters['cnvdelgermline'])
    print('Number of somatic deletion CNVs : ',parameters['cnvdelsomatic'])
    print('Number of somatic aneuploid events : ',parameters['aneuploid'])


    #Get germline variants ---------------------------------------------------------------------------------------------------

    #create empty lists/dictionarys
    germlinevariants=[[],{},{},{},{}]  #variants-0, haplotypes-1, dict of SNV positions-2, dict of CNV breakpoints-3, dict of deleted regions-4

    #get starting chromosome haplotypes
    for chro in gen:
        print(chro)
        germlinevariants[1][chro]=['A','B']
    #add empty keys for each haplotype in each list
        for hap in germlinevariants[1][chro]:
            for l in [2,3,4]:
                germlinevariants[l][chro+hap]=[]

    #get variants
    vartypelist=['cnvrep']*parameters['cnvrepgermline']+['cnvdel']*parameters['cnvdelgermline']+['indel']*indgernum+['snv']*snvgernum
    if len(vartypelist)!=0:
        vartypelist=shufflelist(vartypelist)
    for vartype in vartypelist:
        if vartype == 'cnvrep' or vartype == 'cnvdel':
            germlinevariants=getcnv(gen,germlinevariants,vartype,'germline')
        elif vartype == 'indel':
            germlinevariants=getind(gen,germlinevariants,dbindels,parameters['dbsnpindelproportion'])
        elif vartype == 'snv':
            germlinevariants=getsnv(gen,germlinevariants,dbsnvs,parameters['dbsnpsnvproportion'])


    #Get somatic variants ---------------------------------------------------------------------------------------------------

    #get total distance of all clones from parents
    totdis=0
    for i in clones.keys():
        totdis+=float(clones[i][0])

    #get number of variants per clone
    for clo in clones:
        clones[clo].extend([0,0,0,0,0])    #repcnvs-2, delcnvs-3, indels-4, snvs-5, aneuploids-6
    for i in range(0,parameters['cnvrepgermline']):     #pick random clones to add cnvrep to
        clo=chooseclone(clones,totdis)
        clones[clo][2]+=1
    for i in range(0,parameters['cnvdelgermline']):     #pick random clones to add cnvdel to
        clo=chooseclone(clones,totdis)
        clones[clo][3]+=1
    for i in range(0,indsomnum):     #pick random clones to add indel to
        clo=chooseclone(clones,totdis)
        clones[clo][4]+=1
    for i in range(0,snvsomnum):     #pick random clones to add snv to
        clo=chooseclone(clones,totdis)
        clones[clo][5]+=1
    for i in range(0,parameters['aneuploid']):     #pick random clones to add aneuploid chromosome to
        clo=chooseclone(clones,totdis)
        clones[clo][6]+=1

    #Get somatic variants
    unsortedclones={} #clones whos files haven't been written yet
    sortedclones={}   #clones whos files have been written
    sortedclones['germline']=''
    variants={}
    for clo in clones:
        unsortedclones[clo]=''
    while len(unsortedclones)!=0:
        for clo in list(unsortedclones.keys()):   #write files for clones for which parent clone's files have been written
            if clones[clo][1] in sortedclones:
                variants[clo]=[[],{},{},{},{}] #variants-0, haplotypes-1, dict of SNV positions-2, dict of CNV breakpoints-3, dict of deleted regions-4, dict of indel positions-5
                if clones[clo][1] != 'germline': #if parent clone isn't germline then copy variants from parent clone
                    for a in variants[clones[clo][1]][0]: variants[clo][0].append(a) #copy variants
                    for a in list(variants[clones[clo][1]][1]): variants[clo][1][a]=list(variants[clones[clo][1]][1][a])  #copy chromosome haplotypes
                    for a in list(variants[clones[clo][1]][2]): variants[clo][2][a]=list(variants[clones[clo][1]][2][a])  #copy dict of SNV/indel positions
                    for a in list(variants[clones[clo][1]][3]): variants[clo][3][a]=list(variants[clones[clo][1]][3][a])  #copy dict of CNV breakpoints
                    for a in list(variants[clones[clo][1]][4]): variants[clo][4][a]=list(variants[clones[clo][1]][4][a])  #copy dict of deleted regions
                else:
                    #get variants from germline
                    for a in germlinevariants[0]: variants[clo][0].append(a)
                    for a in list(germlinevariants[1]): variants[clo][1][a]=list(germlinevariants[1][a])
                    for a in list(germlinevariants[2]): variants[clo][2][a]=list(germlinevariants[2][a])
                    for a in list(germlinevariants[3]): variants[clo][3][a]=list(germlinevariants[3][a])
                    for a in list(germlinevariants[4]): variants[clo][4][a]=list(germlinevariants[4][a])
                #get new varaints
                vartypelist=['cnvrep']*clones[clo][2]+['cnvdel']*clones[clo][3]+['indel']*clones[clo][4]+['snv']*clones[clo][5]+['aneu']*clones[clo][6]
                if len(vartypelist)!=0:
                    vartypelist=shufflelist(vartypelist)
                else:
                    vartypelist=[]
                for vartype in vartypelist:

                    if vartype == 'cnvrep' or vartype == 'cnvdel':
                        variants[clo]=getcnv(gen,variants[clo],vartype,'somatic')
                    elif vartype == 'indel':
                        variants[clo]=getind(gen,variants[clo],'','')
                    elif vartype == 'snv':
                        variants[clo]=getsnv(gen,variants[clo],'','')
                    elif vartype == 'aneu':
                        variants[clo]=getaneu(gen,variants[clo])

                del unsortedclones[clo]
                sortedclones[clo]=''


    #Write variant files------------------------------------------------------------------------------------------------------------

    variants['germline']=germlinevariants[:]
    clones['germline']=''
    writevariantfile(parameters['directory'],parameters['prefix'],germlinevariants,'germline')
    for clo in clones:
        writevariantfile(parameters['directory'],parameters['prefix'],variants[clo],clo)
        writeaneuploidfile(parameters['directory'],parameters['prefix'],variants[clo],clo,'somatic',gen)

    with open(parameters['directory'] + '/' + parameters['prefix'] + 'variants_json.txt','w+') as file:
        json.dump(variants, file, indent=1)


# If run as main, run main():
if __name__ == '__main__': main()
