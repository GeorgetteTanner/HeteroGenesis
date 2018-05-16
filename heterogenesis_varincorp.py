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
from signal import signal, SIGPIPE, SIG_DFL
import json
from copy import deepcopy
import os.path

signal(SIGPIPE, SIG_DFL) # Handle broken pipes

version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

def main():

    parser = argparse.ArgumentParser(description="Create random SNVs, indels and CNVs for each subclone in a tumour sample.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(version['__version__']))
    parser.add_argument('-j', '--json', dest='jsonfile', required=True, type=str, help='Json file with parameters')
    parser.add_argument('-c', '--clone', dest='clone', required=True, type=str, help='Clone to be processed')

    args = parser.parse_args()
    clo=args.clone

    with open(args.jsonfile,'r') as file:
        parameters=json.load(file)

    #set default parameter values and give error/warning messages
    if "prefix" not in parameters:
        parameters['prefix'] = ''
    if "reference" not in parameters:
        print('Error: No input genome fasta file provided.')
    if not os.path.exists(parameters['reference'][:-3] + 'i'):
        print('Error: No fai index file for genome.')    
    if "directory" not in parameters:
        print('Warning: No output directory given, using current directory.')
        parameters['directory']='./'
    if "chromosomes" not in parameters:
        print('Warning: No chromosomes given, using "all".')
        parameters['chromosomes']='all'

    #Functions for reading in data----------------------------------------------------------------------------------

    def readinfai(chromosomes,fai,referencefile): #reads in reference and fai files for required chromosomes into dictionaries
        if chromosomes==['all']:
            keepchromos=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
        else:
            keepchromos=chromosomes
        with open(fai,'r') as geno:
                gen = dict([(line.strip().split("\t")[0], float(line.strip().split("\t")[1])) for line in geno])
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

    def readinvars(parameters):    #formats structure parameter into dictionary
        with open(parameters['directory'] + '/' + parameters['prefix'] + 'variants_json.txt','r') as file:
                variants=json.load(file)
        return variants

    #Functions for writing output files--------------------------------------------------------------------------------

    def writeblocksfile(directory,prefix,clo,hapvars,modchros):
        with open(directory + '/' + prefix + clo + 'blocks.txt','w+') as file:
            for chro in hapvars:
                for hap in hapvars[chro]:
                    file.write(clo+'_'+chro+'_'+hap+'\n')
                    for b in modchros[chro][hap].allblocks:
                        file.write(str(b)+'\n')

    def writecnvfile(directory,prefix,clo,combcnvs):
        with open(directory + '/' + prefix + clo + 'cnv.txt','w+') as file:
            for chro in combcnvs:
                for b in combcnvs[chro]:
                    file.write(chro+'\t'+str(b.start)+'\t'+str(int(b.end))+'\t'+str(b.content)+'\n')

    def writevcffile(directory,prefix,clo,combvcfs):
        with open(directory + '/' + prefix + clo + '.vcf','w+') as file:
            for chro in combvcfs:
                for v in combvcfs[chro]:
                    #write: chromosome, position, ref base, alternate base, frequency, total copies, copies per haplotype, copynumber at position
                    file.write(chro+'\t'+str(combvcfs[chro][v][0])+'\t'+str(combvcfs[chro][v][1])+'\t'+str(combvcfs[chro][v][2])+'\t'+str(combvcfs[chro][v][3])+'\t'+str(combvcfs[chro][v][4])+'\t'+str(combvcfs[chro][v][5])+'\t'+str(combvcfs[chro][v][6])+'\n')

    #Functions for generating output file data-------------------------------------------------------------------------------------

    def countvcfs(branches,n):
        for b in branches:
            if b.content=='var':
                n+=1
            else:
                n=countvcfs(b.content,n)
        return n

    def getstart(block):
        return(block.start)

    class BLOCK(object):
        def __init__(self, start, end, content):
            self.start = start
            self.end = end
            self.content=content
        def includes(self, other):  #self completely includes other
            if (other.start >= self.start) and (other.end <= self.end): return True
            return False
        def splitby(self,other):    #other start or end in middle of self
            if (other.end < self.start) or (other.start > self.end): return False   #speeds up code having this check first
            if (other.start > self.start) and (other.start <= self.end): return True
            elif (other.end >= self.start) and (other.end < self.end): return True
            return False
        def __str__(self): return('BLOCK: {}-{}, {}'.format(self.start, self.end, self.content))

    class VCFVAR(object):
        def __init__(self, pos, ref, alt, branches,final,haplo):
            self.pos = pos
            self.ref = ref
            self.alt = alt
            self.branches = branches
            self.final = final
            self.haplo = haplo
        def incnv(self,cnv):
            if (cnv.start <= self.pos) and (cnv.end >= self.pos): return True
            return False
        def __str__(self): return('VCF: {}, {}, {}, {}, {}'.format(self.pos, self.ref, self.alt, self.branches, self.final, self.haplo))

    class CNVBRANCH(object):
        def __init__(self, start, end, content):
            self.start = start
            self.end = end
            self.content = content
        def withincnv(self,cnv):
            if (cnv.start <= self.start) and (cnv.end >= self.end): return True
            return False
        def __str__(self): return('{}, {}, {}'.format(self.start, self.end, self.content))

    class MODCHRO(object):
        def __init__(self,chromosome,allblocks,cnblocks,vcfcounts):
            self.chromosome=chromosome
            self.allblocks=allblocks
            self.cnblocks=cnblocks
            self.vcfcounts=vcfcounts
        def getbasestring(self,ref):
            basestring=[]
            for b in self.allblocks:
                if b.content=='ref':
                    basestring.extend(reference[self.chromosome][int(b.start)-1:int(b.end)])    #b.start and b.end sometimes get .0 on end so need to be converted to int
                else:
                    basestring.extend(b.content)
            return basestring
        def __str__(self): return('MODCHRO from {}: {} allblocks, {} cnblocks'.format(self.chromosome, len(self.allblocks), len(self.cnblocks)))
        def updateblocks(self,var):
            if type(var.content)==str:  #if var contains a string (ie. snv or insertion indel) then just insert into allblocks
                for b in self.allblocks:
                    if (var.start >= b.start) and (var.end <= b.end):   #b.includes(var) #runs quicker having the code here
                        insertblock(self.allblocks,b,var)  #split block in half, remove one base form middle, and insert new block
                        break    #break is important to stop multiple copies of the block having the vairant inserted
            elif type(var.content)==int:   #if var contains an integer (ie. cnv or deletion indel) then update both allblocks and cnblocks
                for blocks in self.allblocks,self.cnblocks:
                    if blocks==self.allblocks or var.end-var.start+1>50:    #prevent indel deletions being written to cnblocks
                        #split blocks at var break points
                        contin=False
                        while contin==False:
                            contin=True
                            for b in blocks:
                                if ((var.start > b.start) and (var.start <= b.end)) or ((var.end >= b.start) and (var.end < b.end)):    #b.splitby(var) #runs quicker having the code here
                                    blocks=splitblocks(blocks,b,var)  #split block at cnv start and/or end position
                                    contin=False   #break the while loop only if run through all bs with no splits
                                    break   #start from beginning of blocks as blocks have now changed
                        #copy blocks withing the var region
                        i=-1
                        copy=[]
                        for b in blocks:
                            i+=1
                            if var.includes(b):
                                copy.append(b)
                            elif copy!=[]: #if there are blocks recorded to copy
                                if var.content!=0:
                                    blocks[i:i]=copy*(var.content - 1)    #insert copied blocks times copy number -1
                                else:   #remove blocks if var is a deletion
                                    for bl in copy:
                                        blocks.remove(bl)
                                break
                            else:
                                pass
            else:
                print("ERROR, 200")
            return self
        def updatevcf(self,var):
            #if var is a vcfvar (ie. snp or indel) then add to modchro vcfcounts list
            if type(var)==VCFVAR:
                self.vcfcounts.append(deepcopy(var))
            #if var is a cnv then adjust numbers of copies of snvs/indels that are located in the region
            elif type(var)==BLOCK:
                for v in self.vcfcounts:    #for each VCFVAR object in MODCHRO.vcfcounts
                    if v.incnv(var):    #if VCFVAR in cnv block
                        v.branches=adjustbranches(v.branches,var)
        def addupfinalvcfs(self):
            for vcfvar in self.vcfcounts:
                vcfvar.final=countvcfs(vcfvar.branches,0)

    def adjustbranches(level,var):  #adjusts record of cnvs overlapping a var in order to calculate how many copies a hap contains. cnvs are recorded in a tree structure with each level coresponding to a cnv position (with the exception of the bottom level which refers to the var position) and each branch on a level refers to a copy.
        if level==[]:   #if copy has been deleted
            return level    #do nothing
        elif level[0].withincnv(var):
            level=[CNVBRANCH(var.start,var.end,deepcopy(level)) for i in range(0,var.content)] #insert new level. need to create new instances of object instead of just copying pointer to the existing instance
            return level
        else:
            level[0].content=adjustbranches(deepcopy(level[0].content),var)  #go to next level down
        return level

    def splitblocks(blocks,spblock,newblock):
        #split block by an overlapping block and return the fragments of the split block
        if newblock.start > spblock.start and newblock.end < spblock.end:
            leftblock=BLOCK(spblock.start,newblock.start-1,spblock.content)
            midblock=BLOCK(newblock.start,newblock.end,spblock.content)
            rightblock=BLOCK(newblock.end+1,spblock.end,spblock.content)
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=[leftblock,midblock,rightblock]
            return blocks
        elif newblock.start > spblock.start:
            leftblock=BLOCK(spblock.start,newblock.start-1,spblock.content)
            rightblock=BLOCK(newblock.start,spblock.end,spblock.content)
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=[leftblock,rightblock]
            return blocks
        elif newblock.end < spblock.end:
            leftblock=BLOCK(spblock.start,newblock.end,spblock.content)
            rightblock=BLOCK(newblock.end+1,spblock.end,spblock.content)
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=[leftblock,rightblock]
            return blocks


    def insertblock(blocks,spblock,newblock):
        #insert a single base block (from snv or insertion indel) into an existing block and return all 3 resulting blocks (or 2 if the blocks start or end at the same position)
        leftblock=BLOCK(spblock.start,newblock.start-1,spblock.content)
        rightblock=BLOCK(newblock.end+1,spblock.end,spblock.content)
        if leftblock.start==newblock.start: #if newblock starts on same base as existing, don't include leftblock
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=deepcopy([newblock,rightblock])
        elif rightblock.end==newblock.end:  #if newblock ends on same base as existing, don't include rightblock
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=deepcopy([leftblock,newblock])
        else:
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=deepcopy([leftblock,newblock,rightblock])
        return blocks

    def combinecnvs(modchro):
        #put all cnv blocks into one list
        allcnvs=[]
        for hap in modchro:
            allcnvs.extend(modchro[hap].cnblocks)
        contin=False
        #split cnv blocks so none overlap
        while contin==False:
            contin2=False
            for cnv in allcnvs:
                contin=False
                if contin2==True:
                    break
                for c in allcnvs:
                    if c.splitby(cnv):
                        allcnvs=splitblocks(allcnvs,c,cnv)  #split c at cnv start and/or end position
                        contin2=True   #causes break from first loop
                        break   #start from beginning of cnvs as cnvs have now changed
                    else:
                        contin=True    #break the while loop only if run through all cnvs with no splits
        #for each cnv block in allcnvs, if its position is not in recorded then add to combined and sum all cnv blocks with same position in allcnvs
        combined=[]
        recorded={}
        for c in allcnvs:
            if c.start not in recorded:
                combined.append(BLOCK(c.start,c.end,sum(s.content for s in allcnvs if s.start==c.start)))
            recorded[c.start]=''
        #sort cnvs
        combined=sorted(combined,key=getstart)
        return combined

    def combinevcfs(modchro,combcnvs):
        #put all vars in dictionaries referenced by position
        allvcfs={}
        for hap in modchro:
            allvcfs[hap]={}
            for v in modchro[hap].vcfcounts:
                allvcfs[hap][v.pos]=v
        #add vars to another dictionary while combining vars with the same position
        combined={}
        for hap in allvcfs:
            for v in allvcfs[hap]:
                if v not in combined:
                    combined[v]=[allvcfs[hap][v].pos,allvcfs[hap][v].ref,allvcfs[hap][v].alt,'round(total/cn,5)','total',[[allvcfs[hap][v].haplo,allvcfs[hap][v].final]],'copynumber']
                else:
                    combined[v][5].append([allvcfs[hap][v].haplo,allvcfs[hap][v].final])
        #Fill in the missing data
        for v in combined:
            #get total number of variant copies copy
            total=sum([i[1] for i in combined[v][5]])
            combined[v][4]=total
            #get copy number
            for cnv in combcnvs:
                if VCFVAR(combined[v][0],'','','','','').incnv(cnv):
                    cn=cnv.content
            combined[v][6]=cn
            #get frequencies
            combined[v][3]=round(total/cn,5)
        return combined

    def writebasestringtofile(parameters,clo,chro,hap,basestring):
        with open(parameters['directory'] + '/' + parameters['prefix'] + clo+chro+hap+'.fasta','w') as file:
            count=0
            file.write('>'+clo+'_'+chro+'_'+hap+'\n')
            for x in basestring:
                if count==80:   #write 80 bases per line
                    file.write('\n')
                    count=0
                file.write(x)
                count+=1
            file.write('\n')

    def createhapvars(clones,gen,variants):
    #create lists of variants by haplotype
        hapvars={}
        for chro in gen:
            hapvars[chro]={}
            for hap in variants[clo][1][chro]:
                hapvars[chro][hap]=[]
                for var in variants[clo][0]:
                    if var[1]==chro:
                        if var[2] in hap:   #if variant haplotype is parent chromosome of or equals the current haplotype, then add var to list
                            hapvars[chro][hap].append(var)
        return hapvars

    def createmodchros(clones,gen,variants):
    #create starting modchros
        modchros={}
        for chro in gen:
            modchros[chro]={}
            for hap in variants[clo][1][chro]:
                modchros[chro][hap]=MODCHRO(chro,[BLOCK(1,gen[chro],'ref')],[BLOCK(1,gen[chro],1)],[]) #starting allblocks is a block of 1-end containing the refernce, and starting cnblocks is a block of 1-end with copy number of 1.
        return modchros


    #Read in variant_dict file and reference genomes------------------------------------------------------------------------------------------------------------

    variants=readinvars(parameters)

    if clo not in variants:
        print(clo + ' not listed in heterogenesis_vargen.py output. Exiting.')
        exit()

    if type(parameters['chromosomes'])==list:
        chromosomes=parameters['chromosomes']
    else:
        chromosomes=[parameters['chromosomes']]

    gen,reference=readinfai(chromosomes,parameters['fai'],parameters['reference'])  #get dictionaries of genome lengths and sequences


    #Generate vcf and cnv output data and write to files-------------------------------------------------------------------------------------------------
    #convert variants to objects and use to update modchros, and then calculate combined vcfs and cnvs

    modchros=createmodchros(clo,gen,variants)
    hapvars=createhapvars(clo,gen,variants)
    combcnvs={}    #dictionary of combined copy numbers for each chromosome
    combvcfs={}    #dictionary of combined vcfs for each chromosome
    for chro in hapvars:
        for hap in hapvars[chro]:
            for var in hapvars[chro][hap]:
                if var[0]=='cnv':
                    modchros[chro][hap].updateblocks(BLOCK(var[3],var[3]+var[4]-1,var[5]))
                    modchros[chro][hap].updatevcf(BLOCK(var[3],var[3]+var[4]-1,var[5]))
                elif var[0]=='indel':
                    if var[7]=='i':
                        modchros[chro][hap].updateblocks(BLOCK(var[3],var[3],var[6]))
                    else:
                        modchros[chro][hap].updateblocks(BLOCK(var[3]+1,var[3]+1+var[4]-1,0))
                    modchros[chro][hap].updatevcf(VCFVAR(var[3],var[5],var[6],[CNVBRANCH(var[3],var[3],'var')],'',hap))
                elif var[0]=='snv':
                    modchros[chro][hap].updateblocks(BLOCK(var[3],var[3],var[5]))
                    modchros[chro][hap].updatevcf(VCFVAR(var[3],var[4],var[5],[CNVBRANCH(var[3],var[3],'var')],'',hap))
                elif var[0]=='aneu':
                    pass    #no need to do anything
            modchros[chro][hap].addupfinalvcfs()
        combcnvs[chro]=combinecnvs(modchros[chro])
        combvcfs[chro]=combinevcfs(modchros[chro],combcnvs[chro])

    #write output files ------------------------------------------------------------------------------------------------------


    #write variant files
    #writeblocksfile(parameters['directory'],parameters['prefix'],clo,hapvars,modchros)
    writecnvfile(parameters['directory'],parameters['prefix'],clo,combcnvs)
    writevcffile(parameters['directory'],parameters['prefix'],clo,combvcfs)

    #Generate genome sequences and write to files
    for chro in hapvars:
        for hap in hapvars[chro]:
            writestring=modchros[chro][hap].getbasestring(gen[chro])
            writebasestringtofile(parameters,clo,chro,hap,writestring)

# If run as main, run main():
if __name__ == '__main__': main()
