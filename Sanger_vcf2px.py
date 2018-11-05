#THIS FILE IS FROM AUGUST 2017, SO MAY BE SLIGHTLY DEPRECIATED
#!/usr/bin/python
'''This python script takes one Sanger imputed vcf files as input, 
removes SNPs with R2<0.8 and MAF>0.01 (options to change), finds the rsID for each SNP, 
and makes output files:

chr*.dosage.txt.gz
samples.txt

For usage, type from command line:
python Sanger_vcf2px.py -h

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
    parser.add_argument('-m', '--maf',
                        help='maf threshold, default 0.01',
                        type=float,
                        default=0.01
		       )
    parser.add_argument('--cpos', '--rsidtocpos',
                        help='outputs c:pos IDs instead of rsIDs, default is rsID',
                        action='store_true'
		       )
    parser.add_argument('-info', '--info',
                        help='IMPUTE2 INFO threshold, default 0.8',
                        type=float,
                        default=0.8
                        )
    parser.add_argument('-dict', '--printdict',
                        help='print out map file made from dictionary',
                        action='store_true'
                        )
    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
chrpath = args.inputdir
c = args.chr
mafthresh = args.maf
r2thresh = args.info
cpos = args.cpos
dict = args.dict
chrfile = chrpath + c + ".vcf.gz"

# get dosage file data
if(os.path.exists(chrpath + 'Sanger_dosages/') == False):
    os.mkdir(chrpath + 'Sanger_dosages/')
if (dict == True):
    mapfile=gzip.open(chrpath + "Sanger_dosages/cpos_rsid_map_chr" + i + ".txt.gz","wb")
    mapfile.write("cpos\trsid\n")
outdosage = gzip.open(chrpath + "Sanger_dosages/chr" + c + ".maf" + str(mafthresh) + ".info" + str(r2thresh) + ".dosage.txt.gz","wb")
for line in gzip.open(chrfile):
    if(line.startswith('##')): #skip lines until field descriptors
        continue
    arr = line.strip().split()
    if(line.startswith('#CHROM')): #only one line should match #CHROM
        ids = arr[9:]
	outdosage.write("chr snp_ID pos ref alt AA_freq " + " ".join(ids) + '\n')
        #split and join ids into FID and IID for PrediXcan
        ids2 = map(lambda x : x.split("_"), ids)
        ids = map(lambda x : ' '.join(x), ids2)
        outsamples = open(chrpath + "Sanger_dosages/samples.txt","w")
        outsamples.write("\n".join(ids))
        outsamples.close()
        continue
    (chr, pos, id, ref, alt, qual, filter, info, format) = arr[0:9]
    if(dict == True):
        mapfile.write(c + ":" + pos + "\t" + id + "\n") 
    if(len(ref) > 1 or len(alt) > 1): #do not output indels, PrediXcan only allows SNPs
        continue
    if(re.search('TYPED',info) == None): #look for 'TYPED' to decide whether to split into 4 or 5
        (refAF, an, ac, infoscore) = info.split(";")
    else:
        (v1, refAF, an, ac, infoscore) = info.split(";")
    r2 = float(infoscore.split("=")[1]) #get IMPUTE2 info score as float
    gt_dosagerow = arr[9:]
    dosagerow = map(lambda x : float(x.split(":")[2]), gt_dosagerow) #lambda function to split each info entry and collect the dosage
    freqalt = round(sum(dosagerow)/(len(dosagerow)*2),4) #calc ALT allele freq (I found that ALT is not always the minor allele)
    if freqalt < 0.5:
        minor = float(freqalt)
    else:
        minor = 1 - float(freqalt)
    if (cpos == True and r2 > r2thresh and minor > mafthresh):
        dosages = ' '.join(map(str,dosagerow))
	id = chr + ':' + pos
        output = 'chr' + chr + ' ' + id + ' ' + pos + ' ' + ref + ' ' + alt + ' ' + str(freqalt) + ' ' + dosages + '\n'
        outdosage.write(output)
    elif(r2 > r2thresh and minor > mafthresh and re.search('rs',id) != None): #only pull SNPs with rsids and default: INFO>0.8, maf>0.01
        dosages = ' '.join(map(str,dosagerow))
        output = 'chr' + chr + ' ' + id + ' ' + pos + ' ' + ref + ' ' + alt + ' ' + str(freqalt) + ' ' + dosages + '\n'
        outdosage.write(output)

outdosage.close()
