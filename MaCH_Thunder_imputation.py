import gzip
import re
import sys
import argparse
import os

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to filter imputed VCF')
    parser.add_argument('-i', '--input',
                        help='input vcf file',
                        required='True'
                        )
    parser.add_argument('-o', '--out',
                        help='output dir',
			                  type=str,
                        required='True'
                        )
    parser.add_argument('-m', '--maf',
                        help='maf threshold, default 0.01',
                        type=float,
                        default=0.01
		       )
    parser.add_argument('--cpos', '--rsidtocpos',
                        help='outputs build:c:pos IDs instead of rsIDs, default is rsID',
                        action='store_true'
		       )
    parser.add_argument('-r2', '--rsq',
                        help='R2 threshold, default 0.8',
                        type=float,
                        default=0.8
                        )
    parser.add_argument('--dict', '--printdict',
                        help='print out map file made from dictionary',
                        action='store_true'
                        )
    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
file = args.input
#c = args.chr
mafthresh = args.maf
r2thresh = args.rsq
cpos = args.cpos
dict = args.dict
#chrfile = chrpath + c + ".vcf.gz"

# get dosage file data
if(os.path.exists(args.out) == False):
    os.mkdir(args.out)
if (dict == True):
    mapfile=gzip.open(args.out +"/cpos_rsid_map.txt.gz","wb")
    mapfile.write("cpos\trsid\n")
outdosage = gzip.open(args.out + "/Mach.maf" + str(mafthresh) + ".R2" + str(r2thresh) + ".dosage.txt.gz","wb")
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
        outsamples = open(args.out + "/samples.txt","w")
        outsamples.write("\n".join(ids))
        outsamples.close()
        continue
    (chr, pos, id, ref, alt, qual, filter, info, format) = arr[0:9]
    if(dict == True):
        mapfile.write(chr + ":" + pos + "\t" + id + "\n") 
    if(len(ref) > 1 or len(alt) > 1): #do not output indels, PrediXcan only allows SNPs
        continue
    try:
        rsqs = re.search('RSQ=0\.[0-9]*',info).group(1) #search for RSQ value
    except AttributeError:
        continue
    r2 = float(rsqs.split("=")[1]) #get IMPUTE2 info score as float
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
    elif(r2 > r2thresh and minor > mafthresh): #only pull SNPs with rsids and default: INFO>0.8, maf>0.01
        dosages = ' '.join(map(str,dosagerow))
        output = 'chr' + chr + ' ' + id + ' ' + pos + ' ' + ref + ' ' + alt + ' ' + str(freqalt) + ' ' + dosages + '\n'
        outdosage.write(output)

outdosage.close()
