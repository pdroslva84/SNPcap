#!/usr/bin/env python

#####################
# gtvalues2plink.py #
#####################
'''
takes genotypes and corresponding GQ, RGQ and DP values (output by GATK's VariantToTable)
and converts to plink .map and .ped, filtering by user-defined values

input has columns (with header):
CHROM, POS, {sample1}.GT, {sample1}.GQ, {sample1}.RGQ, {sample1}.DP, {sample2}...

usage:
$ python gtvalues2plink.py <input.gt.values> [-minGQ -minRGQ -minDP]
'''

import sys
import argparse

# parsing command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument('genotypes', help='sample genotypes with GQ, RGQ and DP vales')
argparser.add_argument('-minGQ', type=int, default=None, help='minimum Genotype Quality (only for non-reference sites')
argparser.add_argument('-minRGQ', type=int, default=None, help='minimum Reference Genotype Quality (only for reference sites')
argparser.add_argument('-minDP', type=int, default=None, help='minimum Depth of Coverage')
args = argparser.parse_args()


if (args.minGQ and not args.minRGQ) or (not args.minGQ and args.minRGQ):
    raise Exception('minGQ and minRGQ have to be defined together!')


output_basename = args.genotypes[:-7]
output_params = '-'.join(['minGQ' + str(args.minGQ),
                          'minRGQ' +  str(args.minRGQ),
                          'minDP' + str(args.minDP)])

with open(args.genotypes, 'r') as inf, \
     open(output_basename + '_' + output_params + '.ped', 'w') as out_ped, \
     open(output_basename + '_' + output_params + '.map', 'w') as out_map:

    sample_list = []
    gt_dict= {}
    pos_list = []

    for line in inf:

        line_list = line.strip().split()

        if line_list[0] == 'CHROM':
            for header in line_list[2::4]:
                sample_list.append(header.split('/')[-1].split('.')[0])
            gt_dict = {sample : [] for sample in sample_list}
            continue

        chrom, pos = line_list[:2]
        snp_name = chrom + '.' + pos
        chrom = chrom.strip('chr')
        pos_list.append([chrom, snp_name, '0', pos])

        values_dict = {}
        for s, sample in enumerate(sample_list):
            gt = line_list[s*4+2]
            gq = line_list[s*4+2+1]
            rgq = line_list[s*4+2+2]
            dp = line_list[s*4+2+3]

            values_dict[sample] = {'GT': gt if gt != './.' else None,
                                   'GQ': int(gq) if gq != 'NA' else None,
                                   'RGQ': int(rgq) if rgq != 'NA' else None,
                                   'DP': int(dp) if dp != 'NA' else None
                                  }

        for sample in sample_list:
            #if (values_dict[sample]['GQ'] > args.minGQ or \
            #   values_dict[sample]['RGQ'] > args.minRGQ) and \
            #   values_dict[sample]['DP'] > args.minDP:

            if (not args.minGQ and not args.minRGQ) and (values_dict[sample]['DP'] >= args.minDP):
                sample_genotype = values_dict[sample]['GT']

            elif (args.minGQ and args.minRGQ) and (not args.minDP):
                if (values_dict[sample]['GQ'] >= args.minGQ) or (values_dict[sample]['RGQ'] >= args.minRGQ):
                    sample_genotype = values_dict[sample]['GT']

            elif (args.minGQ and args.minRGQ) and (args.minDP):
                if (values_dict[sample]['GQ'] >= args.minGQ or \
                    values_dict[sample]['RGQ'] >= args.minRGQ) and \
                    values_dict[sample]['DP'] >= args.minDP:

                    sample_genotype = values_dict[sample]['GT']

            else:
                sample_genotype = None

            if not sample_genotype:
                plink_genotype = '0\t0'
            elif len(sample_genotype) > 3:
                plink_genotype = '0\t0'
            else:
                plink_genotype = '\t'.join(sample_genotype.split('/'))

            gt_dict[sample].append(plink_genotype)


    for sample in sample_list:
        sample_lead = '0\t%s\t0\t0\t0\t0\t' % sample
        sample_genotypes = '\t'.join(gt_dict[sample])
        out_ped.write(sample_lead + sample_genotypes + '\n')

    for snp in pos_list:
        out_map.write('\t'.join(snp) + '\n')
