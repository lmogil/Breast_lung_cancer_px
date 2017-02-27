#!/usr/bin/python
'''This python script takes the UMich imputed vcf files as input, 
removes SNPs with R2<0.8, finds the rsID for each SNP, and makes output files 
for each autosome for PrediXcan:
chr${i}.dosage.txt.gz
samples.txt
dose allele is Allele2, see https://github.com/hakyimlab/PrediXcan/blob/master/Software/HOWTO-beta.md'''
from __future__ import division
import gzip
import re
import sys


chrpath = "/home/lmogil/data/lung_cancer_imp"
#chrpath = sys.argv[1]
c = "22"
#c = sys.argv[2]

print c

chrfile = chrpath+"/chr_"+c+"/chr"+c+".dose.vcf.gz"
print chrfile

##make dictionary: keys->positions values->rsids
snpfile = "/home/lmogil/data/breast_cancer_imp/hrcbychr/chr"+c+"_HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
posdict = {}
for line in gzip.open(snpfile):
    arr = line.strip().split()
    posdict[arr[1]] = arr[2]

# get dosage file data
outdosage = gzip.open("/home/lmogil/data/chr"+c+ ".dosage.txt.gz","wb")
for line in gzip.open(chrfile):
	if(line.startswith('##')):
		continue
	arr = line.strip().split()
		#(CHROM, POS, ID, REF, ALT,QUAL, FILTER, INFO, FORMAT) = arr[0:9]
    	if(line.startswith('#CHROM')): #only one line should match #CHROM
    		ids = arr[9:]
    		ids2 = map(lambda x : x.split("_"), ids)
    		ids = map(lambda x : ' '.join(x), ids2)
    		outsamples = open("samples.txt","w")
    		outsamples.write("\n".join(ids))
    		outsamples.close()
    		continue
        (chr, pos, id, ref, alt, qual, filter, info, format) = arr[0:9]
    	if(re.search('ER2',info) == None): #look for 'ER2' to decide whether to split into 3 or 4
        	(af2, maf, impr2) = info.split(";")
    	else:
        	(af, maf, impr2, imper2) = info.split(";")
        r2 = float(impr2.split("=")[1]) #get r2 value as float
        minor = float(maf.split("=")[1])
        rsid = posdict[pos]
    	if(r2 >= 0.8 and minor > 0.1 and re.search('rs',rsid) != None): #only pull SNPs with rsids and R2>=0.8
        	gt_dosagerow = arr[9:]
        #see http://www.python-course.eu/lambda.php for details
        	dosagerow = map(lambda x : float(x.split(":")[1]), gt_dosagerow) #lambda function to split each info entry and collect the dosage
        	freqalt = round(sum(dosagerow)/(len(dosagerow)*2),4) #calc ALT allele freq (I found that ALT is not always the minor allele)
        	dosages = ' '.join(map(str,dosagerow))
        	output = 'chr' + chr + ' ' + rsid + ' ' + pos + ' ' + ref + ' ' + alt + ' ' + str(freqalt) + ' ' + dosages + '\n'
        	outdosage.write(output)

outdosage.close()
