import sys
import math
import bisect
import numpy as np
from optparse import OptionParser


def CheckLD(Z, ld, ld_thresh, z_min):
   high_Z = [i for i in range(len(Z)) if abs(Z[i])>z_min]
   mismatches= []
   if len(high_Z) > 0:
        for i in range(len(high_Z)):
            for j in range(i+1,len(high_Z)):
                x = high_Z[i]
                y = high_Z[j]
                ld_xy = ld[x,y]
                if abs(ld_xy) > ld_thresh:
                    if Z[x]*Z[y] > 0 and ld_xy<0:
                        mismatches.append((x,y))
                    elif Z[x]*Z[y] < 0 and ld_xy>0:
                        mismatches.append((x,y))
   return mismatches

def Read_Locus(filename, Z_header):
    locus = open(filename)
    header = locus.readline().strip().split()
    z_ind = [i for i in range(len(header)) if header[i]==Z_header]
    z_vector =[]
    if(len(z_ind) == 0):
        print("Supplied Z-score header not found. Please use the -z flag to re-specify.")
        sys.exit()
    all_info =[]
    j = z_ind[0]
    for snps in locus:
        split = snps.strip().split()
        all_info.append(snps)
        z = split[j]
        if split[j] == 'NA':
            z = 0
        z_vector.append(float(z))
    locus.close()
    return [all_info, z_vector]
    

def main():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="filename")
    parser.add_option("-t", "--threshold", dest="threshold", default=".7")
    parser.add_option("-l", "--ld_suffix", dest="ld_suff", default="ld")
    parser.add_option("-m", "--min_z", dest="z_min", default =3)
    parser.add_option("-z", "--zheader", dest="z_header", default ="Zscore")
    parser.add_option("-o", "--outname", dest="outname", default ="LD.errors")
    (options, args) = parser.parse_args()

    fname= options.filename
    thresh = float(options.threshold)
    ld_suff = options.ld_suff
    z_min = float(options.z_min)
    head = options.z_header
    outname=options.outname
    usage = \
    """Usage:
        -i (required) specify input file with listing loci
        -t (optional) specify LD threshold (defualt: 0.7)
        -l (optional) specify suffix for ld file (default = "ld")
        -m (optional) specify minimum z-score to compare (default =3)
        -z (optional) specify z score header in locus file (default = Zscore) 
        -o (optional) specify output file name (default = LD.errors)"""

    if(fname == None):
        sys.exit(usage)
    
    locus_files = open(fname)
    flag=0
    out =open(outname, 'w')
    out.write("Loci with potential errors: \n")

    for files in locus_files:
        print("Checking Locus: " +files)
        [locus_info, locus_z] = Read_Locus(files.strip(), head)
        ld_in = np.loadtxt(files.strip()+"."+ld_suff)
        detected_mismatches = CheckLD(locus_z, ld_in, thresh, z_min)
        if(len(detected_mismatches)>0):
            out.write(files)
        for mismatches in detected_mismatches:
            flag = 1
            out.write("Warning: Detected Potential LD error between: \n" +
            locus_info[mismatches[0]] +locus_info[mismatches[1]] +
            "LD between snps in data: " + str(ld_in[mismatches[0],mismatches[1]])+"\n")
    out.close()    
    if(flag ==1):
        print("""LD errors detected. Please see LD.errors file for details.\n
        Please double check to ensure that reference alleles in LD reference panel match reference alleles in GWAS study""")
    else:
        print("No mismatches detected")
if __name__ == "__main__": main()
                    

    

