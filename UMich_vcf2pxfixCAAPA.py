#!/usr/bin/python
'''This python script takes one UMich imputed vcf files as input, 
removes SNPs with R2<0.8 and MAF>0.01 (options to change), finds the rsID for each SNP, 
and makes output files:

chr*.dosage.txt.gz
samples.txt

For usage, type from command line:
python UMich_vcf2px.py -h

dose allele is Allele2, see https://github.com/hakyimlab/PrediXcan/blob/master/Software/HOWTO-beta.md'''

import gzip
import re
import sys
import argparse
import os

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to filter imputed VCF')
    parser.add_argument('-i', '--inputdir',
                        help='directory containing VCF file',
                        required='True'
                        )
    parser.add_argument('-c', '--chr',
                        help='chromosome',
                        type=str,
                        required='True'
                        )
    parser.add_argument('-r', '--refpop',
			help='reference population, hrc or 1000g',
			type=str,
			required='True'
			)
    parser.add_argument('-m', '--maf',
                        help='maf threshold, default 0.01',
                        type=float,
                        default=0.01)
    parser.add_argument('-r2', '--rsq',
                        help='R2 threshold, default 0.8',
                        type=float,
                        default=0.8
                        )
    parser.add_argument('-o', '--outdir',
                        help='Output directory name',
                        type=str,
                        default="Mich"
                        )
    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
chrpath = args.inputdir
c = args.chr
refpop = args.refpop
mafthresh = args.maf
r2thresh = args.rsq
out_dir = args.outdir
if out_dir.endswith("/"):
    out_dir = out_dir
else:
    out_dir = out_dir + "/"

print(c)

chrfile = chrpath + "chr" + c + ".dose.vcf.gz"

##make dictionary: keys->positions values->rsids
posdict = {}
if(refpop == 'hrc'):
    snpfile = "/home/peter/AA_nonGAIN_SCZ/Imputation/ReferencePanel/chr" + c + "_HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
    for line in gzip.open(snpfile):
        arr = line.strip().split()
        posdict[arr[1]] = arr[2]
elif(refpop == '1000g'):
    snpfile = "/home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr"+c+".legend.gz"
    for line in gzip.open(snpfile):
        if(line.startswith('id')):
            continue
        arr = line.strip().split()
        if(arr[0].count(":") == 3): #exclude CNV
          (rs, pos, a0, a1) = arr[0].split(":")
          cpos = str(c) + ":" + str(pos)
          if(rs.startswith("rs")):
            posdict[cpos] = rs
            #print(rs)
elif(refpop == 'cappa'):
    snpfile = "/home/lauren/ref_impute/all.caapa.sorted.txt"
    for line in open(snpfile):
        if(line.startswith('CHROM')):
            continue
        arr = line.strip().split()
        (chr,pos, rs, a0, a1) = arr[0:5]
        cpos=str(chr)+":"+pos
        posdict[cpos] = rs
else:
    print('need correct refpop:cappa, hrc, or 1000g')

#print(posdict)

# get dosage file data
if(os.path.exists(out_dir) == False):
    os.mkdir(out_dir)
print(c)
outdosage = gzip.open(out_dir + "chr" + c + ".dosage.txt.gz","wb")
for line in gzip.open(chrfile):
    if(line.startswith('##')):
            continue
    arr = line.strip().split()
    #print(line)
    if(line.startswith('#CHROM')): #only one line should match #CHROM
        ids = arr[9:]
        outdosage.write("chr snp_ID pos ref alt " + " ".join(ids) + '\n')
        #split and join ids into FID and IID for PrediXcan
        ids2 = map(lambda x : x.split("_"), ids)
        ids = map(lambda x : ' '.join(x), ids2)
        outsamples = open(out_dir + "samples.txt","w")
        outsamples.write("\n".join(ids))
        outsamples.close()
        continue
    (c, pos, id, ref, alt, qual, filter, info, format) = arr[0:9]
    cpos = str(c) + ":" + str(pos)
    #print(str(info))
    if(bool(re.search('ER2',info)) == True): #look for 'ER2' to decide whether to split into 3 or 4
        (af, maf, impr2, imper2) = info.split(";")
    
    elif(bool(re.search('R2',info)) == True):
        (af, maf, impr2) = info.split(";")
        #print(str(impr2))
    else:
        (af, maf) = info.split(";") #GENOTYPED_ONLY SNPs
        impr2='0=0'
    r2 = float(impr2.split("=")[1]) #get r2 value as float
    minor = float(maf.split("=")[1]) #get maf as float
    #print(cpos)
    if cpos in posdict:
        rsid=posdict[cpos]
        #print(rsid)
    else:
        rsid = '.'
    #print(str(r2) + " " + str(minor) + " " + str(rsid))
    if(r2 > r2thresh and minor > mafthresh and re.search('rs',rsid) != None): #only pull SNPs with rsids and default: R2>0.8, maf>0.01
        gt_dosagerow = arr[9:]
        #see http://www.python-course.eu/lambda.php for details
        dosagerow = map(lambda x : float(x.split(":")[1]), gt_dosagerow) #lambda function to split each info entry and collect the dosage
        freqalt = round(sum(dosagerow)/(len(dosagerow)*2),4) #calc ALT allele freq (I found that ALT is not always the minor allele)
        dosages = ' '.join(map(str,dosagerow))
        output = c + ' ' + rsid + ' ' + pos + ' ' + ref + ' ' + alt + ' ' + str(freqalt) + ' ' + dosages + '\n'
        outdosage.write(output)

outdosage.close()
