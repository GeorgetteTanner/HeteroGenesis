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
import datetime
from sys import stderr, exit

signal(SIGPIPE, SIG_DFL) # Handle broken pipes


version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

def main():
    parser = argparse.ArgumentParser(description="Create random SNVs, indels and CNVs for each subclone in a tumour sample.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(version['__version__']))
    parser.add_argument('-c', '--clones', dest='clonefile', required=True, type=str, help="File with clone proportions in format: 'clone name' \t 'fraction'.")
    parser.add_argument('-d', '--directory', dest='directory', required=True, type=str, help='Directory containing VCF and CNV files.')
    parser.add_argument('-p', '--prefix', dest='prefix', required=True, type=str, help='Prefix of VCF and CNV file names.')
    parser.add_argument('-n', '--name', dest='name', required=True, type=str, help='Output name of tumour sample.')

    args = parser.parse_args()
    
    def warning(msg):
        print('WARNING: {}'.format(msg), file=stderr)
    
    def error(msg, exit_code=1):
        print('ERROR: {}'.format(msg), file=stderr)
        exit(exit_code)

    def getstart(block):
        return(block.start)

    class BLOCK(object):
        def __init__(self, start, end, content, aallele, ballele):
            self.start = start
            self.end = end
            self.content=content
            self.aallele=aallele
            self.ballele=ballele
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
            leftblock=BLOCK(spblock.start,newblock.start-1,spblock.content,spblock.aallele,spblock.ballele)
            midblock=BLOCK(newblock.start,newblock.end,spblock.content,spblock.aallele,spblock.ballele)
            rightblock=BLOCK(newblock.end+1,spblock.end,spblock.content,spblock.aallele,spblock.ballele)
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=[leftblock,midblock,rightblock]
            return blocks
        elif newblock.start > spblock.start:
            leftblock=BLOCK(spblock.start,newblock.start-1,spblock.content,spblock.aallele,spblock.ballele)
            rightblock=BLOCK(newblock.start,spblock.end,spblock.content,spblock.aallele,spblock.ballele)
            blocks[blocks.index(spblock):blocks.index(spblock)+1]=[leftblock,rightblock]
            return blocks
        elif newblock.end < spblock.end:
            leftblock=BLOCK(spblock.start,newblock.end,spblock.content,spblock.aallele,spblock.ballele)
            rightblock=BLOCK(newblock.end+1,spblock.end,spblock.content,spblock.aallele,spblock.ballele)
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
                combined.append(BLOCK(c.start,c.end,sum(s.content for s in cnvs if s.start==c.start),sum(s.aallele for s in cnvs if s.start==c.start),sum(s.ballele for s in cnvs if s.start==c.start)))
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
        warning('Clone proportions add up to ' + str(tot) + ', not 1.')
    #check clones exist
    for clo in clones:
        if not (os.path.exists(args.directory +'/'+ args.prefix + clo + 'cnv.txt') + os.path.exists(args.directory +'/'+ args.prefix + clo + '.vcf')):
            error('Variant profiles for '+clo+' do not exist.')
            
    #Get all cnvs from cnv files
    allcnvs={}
    for clo in clones:
        with open(args.directory +'/'+ args.prefix + clo + 'cnv.txt','r') as file:
            file.readline()
            for line in file:
                cnv=line.strip().split('\t')
                if cnv[0] not in allcnvs:
                    allcnvs[cnv[0]]=[]
                allcnvs[cnv[0]].append(BLOCK(int(cnv[1]),int(cnv[2]),float(cnv[3])*float(clones[clo]),float(cnv[4])*float(clones[clo]),float(cnv[5])*float(clones[clo])))

    #combine all cnvs
    comcnvs={}
    for chromo in allcnvs:
        comcnvs[chromo]=combinecnvs(allcnvs[chromo])

    #write file 
    with open(args.directory + '/'+ args.prefix + args.name + 'cnv.txt','w') as file:
        file.write('Chromosome\tStart\tEnd\tCopy Number\tA Allele\tB Allele\n')
        for chromo in comcnvs:
            for cnv in comcnvs[chromo]:
                file.write(chromo+'\t'+str(cnv.start) +'\t'+str(cnv.end)+'\t'+str(cnv.content)+'\t'+str(cnv.aallele)+'\t'+str(cnv.ballele)+'\n')
    

    #Get all vars from vcf files
    allvars=[]
    for clo in clones:
        with open(args.directory +'/'+ args.prefix + clo + '.vcf','r') as file:
            for line in file:
                if line.startswith('#'): 
                    if line.startswith('##reference=file:'):
                        reference = line[17:].strip('\n')
                    continue
                var=line.split('\t')
                allvars.append([var[0],str(var[1]),var[3],var[4],int(var[9].split(':')[1]),int(var[9].split(':')[4]),float(clones[clo]),str(var[9].split(':')[2])])
                
                
    #combine all vars
    comvars={}
    for var in allvars:
        if var[0]+var[1] in comvars:
            comvars[var[0]+var[1]][4]=comvars[var[0]+var[1]][4]+(float(var[4])*float(var[6]))   #add to current value
        else:
            comvars[var[0]+var[1]]=[var[0],var[1],var[2],var[3],(float(var[4])*float(var[6])),'',var[7][0]]  #multiply number of copies by clone proportion
    
    for var in comvars:
        cn=''
        for cnv in comcnvs[comvars[var][0]]:
            if int(comvars[var][1])>=cnv.start and int(comvars[var][1])<=cnv.end:
                cn=cnv.content
                break
        comvars[var][5]=cn  #copy number 
        
    with open(args.directory +'/'+ args.prefix + args.name + '.vcf','w+') as file:
        file.write('##fileformat=VCFv4.2'+'\n')
        file.write('##fileDate='+str(datetime.datetime.today().strftime('%Y%m%d'))+'\n')
        file.write('##source=heterogenesis_varincorp-'+clo+'\n')
        file.write('##reference=file:'+reference+'\n')
        file.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
        file.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alt allele frequency">\n')
        file.write('##FORMAT=<ID=TC,Number=1,Type=Integer,Description="Total copies of alt allele">\n')
        file.write('##FORMAT=<ID=CN,Number=2,Type=Integer,Description="Copy number at position">\n')
        file.write('##FORMAT=<ID=PH,Number=1,Type=Integer,Description="Phase">\n')
        file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+args.name+'\n')
        for v in comvars:
            #write: chromosome, position, ., ref base, alternate base, ., ., 1, FORMAT,frequency, total copies, copynumber at position
            file.write(comvars[v][0]+'\t'+str(comvars[v][1])+'\t.\t'+str(comvars[v][2])+'\t'+str(comvars[v][3])+'\t.\t.\tNS=1\tAF:TC:CN:PH\t'+str(round(float(comvars[v][4])/float(comvars[v][5]),5))+':'+str(round(comvars[v][4],5))+':'+str(round(comvars[v][5],5))+':'+str(comvars[v][6])+'\n')

# If run as main, run main():
if __name__ == '__main__': main()
