#!/usr/bin/env python

######################
# SNP_concordance.py #
######################

'''
calculates concordance between two PLINK datasets (.map & .ped)

usage:
$ python SNP_concordance.py dataset1_basename dataset2_basename
'''

import sys

plink_base_1 = sys.argv[1]
ped1 = plink_base_1 + '.ped'
map1 = plink_base_1 + '.map'

plink_base_2 = sys.argv[2]
ped2 = plink_base_2 + '.ped'
map2 = plink_base_2 + '.map'

def make_coord_list_from_map(map_file):
    '''
    returns a list of SNPs from a .map file
    to avoid mismatches because of different names, SNPs are given a new designation chromo.position
    '''
    coord_list = []
    with open(map_file, 'r') as inf:
        for line in inf:
            chromo, snp, cM, pos = line.strip().split()
            coord_list.append(chromo + '.' + pos)

    return coord_list

def make_geno_dict_from_ped(ped_file):
    '''
    returns a dict of genotypes from a .ped file
    loci in genotype list will be the same as in .ped (i.e. same as in .map)
    {sample: [(allele1, allele2), (allele1, allele2), ...]}
    '''
    geno_dict = {}
    with open(ped_file, 'r') as inf:
        for line in inf:
            line_list = line.strip().split()
            sample = line_list[1]
            geno_list = [(i, k) for i, k in zip(line_list[6::2], line_list[7::2])]
            geno_dict[sample] = geno_list

    return geno_dict

def compare_genotypes(genotype1, genotype2):
    '''
    compares two genotypes in the form (A1, A2) and (A3, A4)
    returns 1 if concordant, 0 if discordant, 2 if ignore (e.g. missing data present)
    '''
    if ('0' in genotype1) or ('0' in genotype2):
        #print '2', genotype1, genotype2
        return 2
    if sorted(genotype1) == sorted(genotype2):
        #print '1', sorted(genotype1), sorted(genotype2)
        return 1
    else:
        #print '0', genotype1, genotype2
        return 0


# read in SNP coordinates
snps1 = make_coord_list_from_map(map1)
snps2 = make_coord_list_from_map(map2)

# read in genotypes
genotypes1 = make_geno_dict_from_ped(ped1)
genotypes2 = make_geno_dict_from_ped(ped2)

for sample in genotypes1:

    concordances = 0
    discordances = 0
    ignores = 0
    not_found = 0

    sample_genotypes = genotypes1[sample]

    for i, snp in enumerate(snps1):

        if snp == '0.0':
            continue
        geno1 = sample_genotypes[i]

        try:
            idx = snps2.index(snp)
            geno2 = genotypes2[sample][idx]
        except ValueError:
            not_found += 1
            continue

        #print i, snp, geno1, idx, geno2
        comparison = compare_genotypes(geno1, geno2)
        if comparison == 0:
            discordances += 1
        elif comparison == 1:
            concordances += 1
        elif comparison == 2:
            ignores += 1
        else:
            print('something went wrong!')

    concordance_rate = float(concordances) / (float(concordances) +  float(discordances))
    total_positions = concordances + discordances + ignores
    print_out = 'sample: {0}, total # positions: {1}, concordance rate: {2}, # positions ignored: {3}, # SNPs not found in 2nd dataset: {4}'.format(sample, total_positions, concordance_rate, ignores, not_found)
    print(print_out)
