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
import os.path

signal(SIGPIPE, SIG_DFL) # Handle broken pipes


version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

def main():
    parser = argparse.ArgumentParser(description="Create random SNVs, indels and CNVs for each subclone in a tumour sample.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(version['__version__']))
    parser.add_argument('-c', '--clones', dest='clonefile', required=True, type=str, help="File with clone proportions in format: 'clone name' \t 'fraction'.")
    parser.add_argument('-d', '--directory', dest='directory', required=True, type=str, help='Directory containing VCF and CNV files.')
    parser.add_argument('-p', '--prefix', dest='prefix', required=True, type=str, help='Prefix of VCF and CNV file names.')


    args = parser.parse_args()

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
            if (other.start > self.start) and (other.start <= self.end): return True
            elif (other.end >= self.start) and (other.end < self.end): return True
            return False
        def __str__(self): return('BLOCK: {}-{}, {}'.format(self.start, self.end, self.content))


    def splitblocks(blocks,spblock,newblock):
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

    def combinecnvs(cnvs):
        contin=False
        while contin==False:
            contin2=False
            for cnv in cnvs:
                contin=False
                if contin2==True:
                    break
                for c in cnvs:
                    if c.splitby(cnv):
                        cnvs=splitblocks(cnvs,c,cnv)  #split c at cnv start and/or end position
                        contin2=True   #causes break from first loop
                        break   #start from beginning of cnvs as cnvs have now changed
                    else:
                        contin=True    #break the while loop only if run through all cnvs with no splits
        combined=[]
        recorded={}
        for c in cnvs:
            if c.start not in recorded:
                combined.append(BLOCK(c.start,c.end,sum(s.content for s in cnvs if s.start==c.start)))
            recorded[c.start]=''
        combined=sorted(combined,key=getstart)
        return combined


    #Read in clone proportions
    clones={}
    with open(args.clonefile,'r') as file:
        for line in file:
            c,p=line.strip().split('\t')
            clones[c]=p
    #Check proportions add up to 1
    tot=sum([float(clones[x]) for x in clones])
    if round(tot,5) != 1:
        print('Warning: Clone proportions add up to ' + str(tot) + ', not 1.')


    #Get all cnvs from cnv files
    allcnvs={}
    for clo in clones:
        with open(args.directory + args.prefix + clo + 'cnv.txt','r') as file:
            for line in file:
                cnv=line.strip().split('\t')
                if cnv[0] not in allcnvs:
                    allcnvs[cnv[0]]=[]
                allcnvs[cnv[0]].append(BLOCK(int(cnv[1]),int(cnv[2]),float(cnv[3])*float(clones[clo])))

    #combine all cnvs
    comcnvs={}
    for chromo in allcnvs:
        comcnvs[chromo]=combinecnvs(allcnvs[chromo])

    #write file 
    with open(args.directory + args.prefix + 'tumourcnv.txt','w') as file:
        for chromo in comcnvs:
            for cnv in comcnvs[chromo]:
                file.write(chromo+'\t'+str(cnv.start) +'\t'+str(cnv.end)+'\t'+str(cnv.content)+'\n')
    

    #Get all vars from vcf files
    allvars=[]
    for clo in clones:
        with open(args.directory + args.prefix + clo + '.vcf','r') as file:
            for line in file:
                var=line.split('\t')
                allvars.append([var[0],str(var[1]),str(var[2]),str(var[3]),int(var[5]),float(clones[clo])]) #(var[5])-number of copies, (clones[clo])-clone proportion 
                
                
    #combine all vars
    comvars={}
    for var in allvars:
        if var[0]+var[1] in comvars:
            comvars[var[0]+var[1]][4]=comvars[var[0]+var[1]][4]+(float(var[4])*float(var[5]))   #add to current value
        else:
            comvars[var[0]+var[1]]=[var[0],var[1],var[2],var[3],(float(var[4])*float(var[5]))]  #multiply number of copies by clone proportion
    for var in comvars:
        for cnv in comcnvs[comvars[var][0]]:
            if int(comvars[var][1])>=cnv.start and int(comvars[var][1])<=cnv.end:
                cn=cnv.content
                break
        comvars[var][4]=round(float(comvars[var][4])/float(cn),5)  #divide total number of copies by overall copy number to get overall VAF   
    with open(args.directory + args.prefix + 'tumour.vcf','w+') as file:
        for var in comvars:
            file.write(comvars[var][0]+'\t'+str(comvars[var][1]) +'\t.\t.\t'+comvars[var][2]+'\t'+str(comvars[var][3])+'\t.\t.\t'+str(round(comvars[var][4],5))+'\n')

# If run as main, run main():
if __name__ == '__main__': main()
