import re

from tqdm.notebook import tqdm

import pandas as pd

# Load the reference GFF file
gff = pd.read_csv('./data/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff', names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'], sep='\t', header=None)
# Load the WHO catalogue
catalogue_1 = pd.read_excel('./data/WHO-UCN-GTB-PCI-2021.7-eng 2.xlsx', sheet_name='Mutation_catalogue', header=[0,1]).set_index([('variant (common_name)', 'Unnamed: 2_level_1')])
# Load the reference genome to impute missing data from deletions
h37rv = ''
f = open('./data/h37rv_reference/h37rv.fasta', 'r')
f.readline()
for line in f.readlines():
    h37rv += line.replace('\n', '')

a_map_1 = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    '*': '!',
}

aa_map_2 = {}
for i in aa_map_1.keys():
    aa_map_2[aa_map_1[i]] = i

nucleotide_complements = {
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'A': 'T',
}

# Setup the regular expressions
re_c = re.compile('^(\w+)_([actg])(-*\d+)([actg])$')
re_p = re.compile('^(\w+)_([A-Z])(\d+)([A-Z!])$')
re_d = re.compile('^(\w+)_(-*\d+)_del_(\d+)_([actg]+)_([actg]+)$')
re_i = re.compile('^(\w+)_(-*\d+)_ins_(\d+)_([actg]+)_([actg]+)$')

def func(row):
    return row.str.lower()

def process_variant(variant):
    '''Translates variants in the WHO catalogue format to HGVS'''
    m = re_c.match(variant)

    if m:
        if gff_dict[m[1]]['type'] == 'rRNA':
            v_type = 'n'
            ref = m[2].upper()
            alt = m[4].upper()
        else:
            v_type = 'c'
            ref = m[2].upper()
            alt = m[4].upper()
        return (m[1], v_type, '{}.{}{}>{}'.format(v_type, m[3], ref, alt), False, None)

    m = re_p.match(variant)
    if m:
        return (m[1], 'p', 'p.{}{}{}'.format(aa_map_2[m[2].upper()], m[3], aa_map_2[m[4].upper()]), False, None)

    m = re_d.match(variant)
    if m:
        if int(m[3]) != len(m[4]) - len(m[5]):
            return (None, None, None, True, 'length mismatch')

        starts = [pos for pos in range(1, len(m[4]) + 1 - int(m[3])) if m[4][:pos]+m[4][pos+int(m[3]):] == m[5]]
        if not starts:
            return (None, None, None, True, 'invalid indel')
        if not gff_dict[m[1]]['strand']:
            hgvs = []
            for start in starts:
                if int(m[3]) == 1:
                    hgvs.append('c.{}del'.format(int(m[2])+start))
                else:
                    hgvs.append('c.{}_{}del'.format(int(m[2])+start, int(m[2])+start-1+int(m[3])))
            return (m[1], 'c', '|'.join(hgvs), False, None)
        else:
            hgvs = []
            for start in starts:
                if int(m[3]) == 1:
                    hgvs.append('c.{}del'.format(int(m[2]) - start - int(m[3]) + 1))
                else:
                    v = 'c.{}_{}del'.format(int(m[2]) - start - int(m[3]) + 1, int(m[2]) - start)
                    hgvs.append(v)
            return (m[1], 'c', '|'.join(hgvs), False, None)

    m = re_i.match(variant)
    if m:
        if int(m[3]) != len(m[5]) - len(m[4]):
            return (None, None, None, True, 'length mismatch')
        starts = [pos for pos in range(1, len(m[4]) + 1) if m[4][:pos]+m[5][pos:pos+int(m[3])]+m[4][pos:] == m[5]]
        if not starts:
            return (None, None, None, True, 'invalid indel')
        if not gff_dict[m[1]]['strand']:
            hgvs = []
            for start in starts:
                hgvs.append('c.{}_{}ins{}'.format(int(m[2])+start-1, int(m[2])+start, ''.join([i.upper() for i in m[5][start:start+int(m[3])]])))
            return (m[1], 'c', '|'.join(hgvs), False, None)
        else:
            hgvs = []
            for start in starts:
                v = 'c.{}_{}ins{}'.format(int(m[2])-start, int(m[2]) - start+1, ''.join([nucleotide_complements[i.upper()] for i in m[5][start:start+int(m[3])][::-1]]))
                hgvs.append(v)
            return (m[1], 'c', '|'.join(hgvs), False, None)

    print(variant)
    return (None, None, None, True, 'does not match indel or variant')

# Get the gene information from the GFF file
info_regex = re.compile('Locus=(.*);Name=(.*);Function=')
gff[['locus_tag', 'name']] = gff.attributes.apply(lambda info: pd.Series(info_regex.match(info).groups([1,2])))

gff_dict = {}
for _, row in gff.iterrows():
    gene = {}
    gene['seqid'] = row.seqid
    gene['source'] = row.source
    gene['type'] = row.type
    gene['start'] = row.start
    gene['end'] = row.end
    gene['score'] = row.score
    gene['strand'] = 0 if row.strand == '+' else 1
    gene['attributes'] = row.attributes
    gene['locus_tag'] = row.locus_tag
    gene['name'] = row['name']
    gff_dict[row['locus_tag']] = gene
    gff_dict[row['name']] = gene


# Prepare the WHO catalogue dataframe
classified = []
v = re.compile('^(.*) \((.*)\)')
for var, row in catalogue_1[catalogue_1[('FINAL CONFIDENCE GRADING', 'Unnamed: 51_level_1')].apply(lambda conf: conf != 'combo')].iterrows():
    drug = row[('drug', 'Unnamed: 0_level_1')]
    m = v.match(var)
    if m:
        # Include all variants listed
        variants = [m[1]] + [i.strip() for i in m[2].split(',')]
        # Include only first variant (the one on which the analysis was performed)
        variants = [m[1]]
    else:
        variants = [var]
    category = row[('FINAL CONFIDENCE GRADING', 'Unnamed: 51_level_1')].split(')')[0]
    genome_pos = '{:.0f}'.format(row[('Genome position', 'Unnamed: 3_level_1')])
    for variant in variants:
        classified.append([variant, drug, category, genome_pos, var])
classified = pd.DataFrame(classified, columns=['variant', 'drug', 'classification', 'genome_position', 'who_original'])

# Convert the variants to HGVS format
for idx, row in tqdm(classified.iterrows(), total=classified.shape[0]):
    x = process_variant(row.variant)
    classified.loc[idx, 'gene'] = x[0]
    classified.loc[idx, 'type'] = x[1]
    classified.loc[idx, 'hgvs'] = x[2]
    classified.loc[idx, 'fail'] = x[3]
    classified.loc[idx, 'fail_reason'] = x[4]

# Impute missing data for deletions
lenght_mismatch = classified[classified.fail_reason == 'length mismatch'].sort_values(by='variant', key=func)

for idx, row in tqdm(lenght_mismatch.iterrows(), total=lenght_mismatch.shape[0]):

    m = re_d.match(row.variant)
    if m:
        if not gff_dict[m[1]]['strand']:
            indexing_correction = -1 if int(m[2]) < 0 else -2 # correct for 0 based python indexing (-1 if promotor, -2 if within gene)
            start = gff_dict[m[1]]['start'] + int(m[2]) + indexing_correction
            end = start + int(m[3]) + len(m[5]) # add the lenght of the alt allele to account for the bases not part of the indel
            complete_variant = '{}_{}_del_{}_{}_{}'.format(m[1], m[2], m[3], h37rv[start:end].lower(), m[5])
            classified.loc[idx, 'complete_variant'] = complete_variant
            classified.loc[idx, 'complete_variant_fail'] = False
        else:
            indexing_correction = -1 if int(m[2]) < 0 else 0 # correct for 0 based python indexing (-1 if promotor, 0 if within gene)
            start = gff_dict[m[1]]['end'] - int(m[2]) + indexing_correction # subtract m[2] instead of adding as this is the opposite strand
            end = start + int(m[3]) + len(m[5]) # add the lenght of the alt allele to account for the bases not part of the indel
            complete_variant = '{}_{}_del_{}_{}_{}'.format(m[1], m[2], m[3], h37rv[start:end].lower(), m[5])
            classified.loc[idx, 'complete_variant'] = complete_variant
            classified.loc[idx, 'complete_variant_fail'] = False
        continue

    m = re_i.match(row.variant)
    if m:
        classified.loc[idx, 'complete_variant_fail'] = True
        classified.loc[idx, 'complete_variant_fail_reason'] = 'Not assuming for insertions'
        if not gff_dict[m[1]]['strand']:
            pass
        else:
            pass
        continue

# Convert imputed deletions to HGVS format
for idx, row in tqdm(classified[classified.complete_variant_fail == False].iterrows(), total=classified[classified.complete_variant_fail == False].shape[0]):
    if row.complete_variant_fail:
        continue
    x = process_variant(row.complete_variant)
    classified.loc[idx, 'gene'] = x[0]
    classified.loc[idx, 'type'] = x[1]
    classified.loc[idx, 'hgvs'] = x[2]
    classified.loc[idx, 'fail'] = x[3]
    classified.loc[idx, 'fail_reason'] = x[4]

# Write results to csv file
classified.to_csv('/Users/lennertverboven/Desktop/who_catalogue.csv', index=False)
