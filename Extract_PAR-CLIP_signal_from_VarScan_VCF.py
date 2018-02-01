# This script extracts the coverage of PAR-CLIP conversion events from the VarScan output VCF files;

import os, sys
from operator import itemgetter

vcf_filenames = sys.argv[1:]
results = []
for vcf_filename in vcf_filenames:
    if os.path.isfile(vcf_filename) and vcf_filename.endswith('.vcf'):
        print(vcf_filename)
        with open(vcf_filename, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                chrom, start, nucl, subst, info = fields[0], int(fields[1]), fields[3], fields[4], fields[8].split(':')
                if nucl == 'T' and subst == 'C':
                    strand = 'fw'
                elif nucl == 'A' and subst == 'G':
                    strand = 'rev'
                else:
                    continue
                end = int(start) + 1
                samples = [sample.split(':') for sample in fields[9:]]
                values = []
                for sample in samples:
                    d  ={}
                    for i in range(len(sample)):
                        d[info[i]] = sample[i]
                    if 'AD' in d:
                        value = int(d['AD'])
                    else:
                        value = 0
                    if strand == 'fw':
                        values.append(value)
                    elif strand == 'rev':
                        values.append(-value)
                results.append([chrom, start, end] + values)
results = sorted(results, key = itemgetter(0, 1, 2))
print(len(results), 'unique conversion positions in total;')

nsamples = len(results[0][3:])
for n in range(nsamples):
    sample_name = 's' + '{:02d}'.format(n+1)
    output_filename = 'PAR-CLIP_signal_' + sample_name + '.bg'
    with open(output_filename, 'w') as output_file:
        output_file.write('track type=bedGraph name=' + sample_name + ' description=' + sample_name + ' color=0,100,200 altColor=200,100,0\n')
        count = 0
        for result in results:
            value = result[n+3]
            if value != 0:
                count += 1
                out = result[:3] + [value]
                output_file.write('\t'.join([str(f) for f in out]) + '\n')
        print(sample_name + ':', count, 'positions;')
print('Done!')