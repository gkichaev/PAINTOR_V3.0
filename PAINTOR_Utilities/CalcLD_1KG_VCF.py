import sys
import numpy as np
import gzip
from optparse import OptionParser



def Read_Locus(file_name, effect_allele, alt_allele, position):
    """Read fine mapping file and return a matrix with [positions, effect allele, alternate allele]"""
    file_stream = open(file_name, 'r')
    header = file_stream.readline().strip().split()
    try:
        effect_index = header.index(effect_allele)
        alt_index = header.index(alt_allele)
        pos_index = header.index(position)
    except(ValueError):
        print("Error: Header is mis-specified. Please check and try-again.")
    all_positions = []
    all_data = []
    all_data.append(header)
    counter = 0
    for lines in file_stream:
        counter += 1
        temp = lines.strip().split()
        #drop ambiguous snps from input data and warn!
        A1_allele = temp[effect_index].upper()
        A0_allele = temp[alt_index].upper()
        if (A1_allele == "A" and A0_allele == "T") or (A1_allele == "T" and A0_allele == "A") or \
         (A1_allele == "G" and A0_allele == "C") or (A1_allele == "C" and A0_allele == "G"):
            print("Warning! Ambiguous SNP (AT or GC) in input locus. Dropping row " + str(counter) + " from data:")
            print(temp)
        else:
            all_positions.append(int(temp[pos_index]))
            all_data.append(temp)
    file_stream.close()
    return [all_positions, all_data]


def Extract_Pop_Ids(file_name, population):
    """Get population ids from mapping file"""
    file_stream = open(file_name, 'r')
    header = file_stream.readline().strip().split()
    all_ids = []
    for lines in file_stream:
        temp = lines.strip().split()
        if (temp[2] == population):
            all_ids.append(temp[0])
    file_stream.close()
    return all_ids

def Filter_VCF_rows(file_name, pos):
    """Extract all the rows in the VCF file that have matching positions in the locus file + header"""
    file_stream = gzip.open(file_name, 'r')
    max_position = max(pos)
    out_rows = []
    pos_set = set()
    [pos_set.add(int(i)) for i in pos]

    pos_present = dict((i, False) for i in pos)

    for lines in file_stream:
        if(lines[0:2] != b"##"):
            if(lines[0:6] == b"#CHROM"):
                out_rows.append(lines.strip().split())
            else:
                line_split = lines.strip().split()
                curr_pos = int(line_split[1])
                if(curr_pos in pos_set):
                    out_rows.append(line_split)
                    pos_present[curr_pos] = True
                elif(curr_pos > max_position):
                    break
        else:
            continue
    file_stream.close()
    return [out_rows, pos_present]

def Extract_Pop_Haps(vcf_rows, pop_ids):
    """Extract (continental) population haplotypes and return a numpy matrix with haplotypes"""
    header = vcf_rows[0]
    header_str = [item.decode("utf-8") for item in header]
    extract_index = [header_str.index(ids) for ids in pop_ids]
    pop_haps =[]
    for rows in vcf_rows[1::]:
        hap_list = [rows[i].decode("utf-8").split("|") for i in extract_index]
        hap = [int(item) for sublist in hap_list for item in sublist]
        pop_haps.append(hap)
    return np.matrix(pop_haps)





def Match_Ref_Panel_Alleles(vcf_rows, effect_allele, alt_allele, position):
    """Determine if reference panel alleles match input data alleles.
        If alleles in data do not match what was extracted from ref panel output the reference panel row is deleted"""
    header = vcf_rows[0]
    A1_index = header.index(b"ALT")
    A0_index = header.index(b"REF")
    pos_index = header.index(b"POS")
    flip_index = []
    out_vcf = []
    out_vcf.append(header)

    for rows in vcf_rows[1::]:
        vcf_pos = int(rows[pos_index])
        vcf_A1 = rows[A1_index].decode("utf-8")
        vcf_A0 = rows[A0_index].decode("utf-8")
        curr_snp = position.index(vcf_pos)
        # check if match
        if vcf_A1 == effect_allele[curr_snp] and vcf_A0 == alt_allele[curr_snp]:
            flip_index.append(0)
            out_vcf.append(rows)
        elif vcf_A1 == alt_allele[curr_snp] and vcf_A0 == effect_allele[curr_snp]:
            flip_index.append(1)
            out_vcf.append(rows)
        else:
            print("Warning! Found alleles " + str(vcf_A1) + " and " + str(vcf_A0) + " in refernce panel\n")
            print("Expecting alleles " + str(effect_allele[curr_snp]) + " and " + str(alt_allele[curr_snp]) + "\n")
    return [out_vcf, flip_index]





def Match_Ref_Panel_Alleles(vcf_rows, input_data, pos_header, effect_header, alt_header, Z_header):


    """Determine if reference panel alleles match input data alleles.
        If alleles in data do not match what was extracted from ref panel output the reference panel row is deleted"""
    header_vcf = vcf_rows[0]
    A1_index = header_vcf.index(b"ALT")
    A0_index = header_vcf.index(b"REF")
    pos_index_vcf = header_vcf.index(b"POS")
    out_vcf = []
    out_vcf.append(header_vcf)

    final_locus = []
    header_input = input_data[0]
    final_locus.append(header_input)

    pos_index = header_input.index(pos_header)
    effect_index = header_input.index(effect_header)
    alt_index = header_input.index(alt_header)
    z_index = header_input.index(Z_header)

    
    position = [int(pos[pos_index]) for pos in input_data[1::]]
    effect_allele = [effect[effect_index] for effect in input_data[1::]]
    alt_allele = [alt[alt_index] for alt in input_data[1::]]


    for rows in vcf_rows[1::]:
        vcf_pos = int(rows[pos_index_vcf])
        vcf_A1 = rows[A1_index].decode("utf-8")
        vcf_A0 = rows[A0_index].decode("utf-8")
        curr_snp = position.index(vcf_pos)
        data_A1 = effect_allele[curr_snp].upper()
        data_A0 = alt_allele[curr_snp].upper()
        # check if reference panel matches input data. drop SNPs that do not have the same alleles in ref panel.
        if vcf_A1 == data_A1 and vcf_A0 == data_A0:
            out_vcf.append(rows)
            final_locus.append(input_data[curr_snp+1])
        elif vcf_A1 == data_A0 and vcf_A0 == data_A1:
            out_vcf.append(rows)
            flip_line = input_data[curr_snp+1]
            zflip = -1*float(flip_line[z_index])
            flip_line[z_index] = str(zflip)
            final_locus.append(flip_line)
        else:
            print("Warning! Found alleles " + str(vcf_A1) + " and " + str(vcf_A0) + " in refernce panel\n")
            print("Expecting alleles " + str(effect_allele[curr_snp]) + " and " + str(alt_allele[curr_snp]) + "\n")
    return [out_vcf, final_locus]


def Write_Output(outname, final_data, computed_ld, drop_mono):
    locus_name = outname+ ".processed"
    locus_out = open(locus_name, 'w')
    ld_name = outname + ".ld"
    if(drop_mono == True):
        poly_index = [i for i in range(len(computed_ld)) if np.isnan(computed_ld[i,i]) != True]
        ld_filt = computed_ld[poly_index].T[poly_index]
        np.savetxt(ld_name,ld_filt, fmt='%1.4e')
        locus_out.write(" ".join(final_data[0]) + "\n")
        for i in poly_index:
            locus_out.write(" ".join(final_data[i+1])+"\n")
    else:
        np.savetxt(ld_name, computed_ld, fmt='%1.4e')
        for lines in final_data:
            locus_out.write(" ".join(lines)+"\n")
    locus_out.close()

def Filter_Missing_SNPs_ref(input_data, present_snps, pos_header):
    filtered_data = []
    header = input_data[0]
    pos_index = header.index(pos_header)
    filtered_data.append(header)
    for snps in input_data[1::]:
        snp_position = int(snps[pos_index])
        if(present_snps[snp_position]):
            filtered_data.append(snps)
        else:
            print("Warning! SNP not found in reference panel :")
            print(" ".join(snps) + "\n")

    return filtered_data

def main():

    ##defaults
    drop_mono=True

    ##
    parser = OptionParser()
    parser.add_option("-l", "--locus", dest="locus")
    parser.add_option("-r", "--reference", dest="reference")
    parser.add_option("-e", "--effect_allele", dest="effect_allele")
    parser.add_option("-a", "--alt_allele", dest="alt_allele")
    parser.add_option("-p", "--population", dest="population")
    parser.add_option("-o", "--out_name", dest="out_name")
    parser.add_option("-i", "--position", dest="position")
    parser.add_option("-m", "--map_file", dest="map_file")
    parser.add_option("-z", "--Zhead", dest="Zhead")
    parser.add_option("-d", "--drop_mono", dest="drop_mono")

    (options, args) = parser.parse_args()

    locus_name = options.locus
    reference =  options.reference
    effect_allele = options.effect_allele
    alt_allele = options.alt_allele
    population = options.population
    out_name = options.out_name
    position = options.position
    map_file = options.map_file
    Zhead = options.Zhead
    drop_mono = options.drop_mono
    usage = \
    """ Need the following flags specified (*)
        Usage:
        --locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position) *
        --reference [-r]  specify reference VCF file corresponding to chromosome of locus *
        --map_file [-m] specify reference map file that maps population ids to individuals *
        --position [-i] specify the name of the field in header corresponding to the SNP positions *
        --effect_allele [-e] specify the name of the field in header corresponding to effect allele (i.e the allele encoded as 1) *
        --alt_allele [-a] specify the name of the field in header corresponding to alternate allele (i.e the allele encoded as 0) *
        --population [-p] specify name of continental population {AFR, AMR, EAS, EUR, SAS} to compute LD with *
        --out_name [-o] specify the stem of the output files *
        --Zhead [-z] specify the name of Zscore field in header *
        --drop_mono [-d] should snps found to be monomorphic in reference panel be dropped? [default: True]
        """

    if(locus_name == None or reference == None or effect_allele == None or alt_allele == None or population == None or out_name == None or position == None):
        sys.exit(usage)

    [all_positions, all_data] = Read_Locus(locus_name,effect_allele,alt_allele,position)
    [vcf_rows, found_postions] = Filter_VCF_rows(reference,all_positions)
    pop_ids = Extract_Pop_Ids(map_file,population)
    filtered_data = Filter_Missing_SNPs_ref(all_data, found_postions, position)
    [final_reference, flipped_data] = Match_Ref_Panel_Alleles(vcf_rows, filtered_data,   position, effect_allele, alt_allele, Zhead)
    pop_haps = Extract_Pop_Haps(final_reference, pop_ids)
    ld_mat = np.corrcoef(pop_haps)
    Write_Output(out_name, flipped_data, ld_mat, drop_mono)


if __name__ == "__main__": main()
