import os
import sys
from Bio import SeqIO
from collections import defaultdict

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
    product = 'NA'
    bgc_length = 'NA'
    complete_status = 'Not Near Contig Edge'
    try:
        with open(bgc_gbk) as obg:
            for line in obg:
                line = line.strip()
                if '/contig_edge="True"' in line:
                    complete_status = 'Near Contig Edge' 
        products = set([])
        with open(bgc_gbk) as obg:
            for rec in SeqIO.parse(obg, 'genbank'):
                for feat in rec.features:
                    if feat.type == 'protocluster':
                        try:
                            products.add(feat.qualifiers.get('product')[0])
                        except:
                            pass
                bgc_length = str(len(rec.seq))

        product_string = ' '.join(products) + ' '
        if len(products) == 1:
            product = list(products)[0]
        elif len(products) == 2 and 'NRPS-like' in products and 'NRPS' in products:
            product = 'NRPS'
        elif len(products) >= 2:
            if 'NRPS ' in product_string and 'PKS' in product_string: 
                product = 'multi-type (with NRPS & PKS)'
            elif 'NRPS ' in product_string:
                product = 'multi-type (with NRPS)'
            elif 'PKS' in product_string:
                product = 'multi-type (with PKS)'
            elif 'NI-siderophore' in products and 'NRPS-metallophore' in products:
                product = 'metallophore'
            else:
                 product = 'multi-type'
            
    except:
        sys.stderr.write('Issues parsing BGC Genbank %s\n' % bgc_gbk)
        raise RuntimeError()
    
    return([product, '|'.join(products), complete_status, bgc_length])

as_results_dir = 'antiSMASH_Results/'
meta_data_file = 'MetaData.txt'
gcf_results = 'BiG-SCAPE_Results/network_files/2023-06-30_15-38-05_hybrids_glocal/mix/mix_clustering_c0.30.tsv'
dereplicated_file = 'redundant.txt'

gcf_mapping = {}
with open(gcf_results) as ogrf:
    for line in ogrf:
        line = line.strip()
        ls = line.split('\t')
        gcf_mapping[ls[0]] = ls[1]

body_site = {}
body_site_type = {}
with open(meta_data_file) as omdf:
    for i, line in enumerate(omdf):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        body_site[ls[0]] = ls[-1].split('(')[0].strip()
        body_site_type[ls[0]] = ls[-1].split('(')[1].strip().split(')')[0]

redundant = set([])
with open(dereplicated_file) as ordf:
    for line in ordf:
        line = line.strip()
        redundant.add(line)

print('\t'.join(['sample_id', 'contig_id', 'bgc_id', 'gcf_id', 'product', 'products_list', 'near_contig_edge_status', 'bgc_length_bp', 'pig_body_site', 'pig_body_site_type']))
for s in os.listdir(as_results_dir):
    if s in redundant: continue
    samp_dir = as_results_dir + s + '/'
    for f in os.listdir(samp_dir):
        if '.region' in f and f.endswith('.gbk'):
            bgc_id = f.split('.gbk')[0]
            sample_id = bgc_id.split('.')[0]
            contig_id = '.'.join(bgc_id.split('.')[:2])
            bgc_gbk = samp_dir + f
            product, products, complete_status, bgc_length = parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk)
            print('\t'.join([sample_id, contig_id, bgc_id, gcf_mapping[bgc_id], product, products, complete_status, bgc_length, body_site[sample_id], body_site_type[sample_id]]))
