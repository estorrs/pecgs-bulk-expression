"""
from https://github.com/ding-lab/HTAN_bulkRNA_expression/blob/master/scripts/gen_fpkm.py
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import mygene


def readcount_to_fpkm(gene_info_df, rc_df):
    # FPKM and FPKM-UQ conversion formula from GDC
    # FPKM: <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm>
    # FPKM-UQ: <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#upper-quartile-fpkm>
    protein_coding_gene_info_df = gene_info_df.query("gene_type == 'protein_coding'")
    gene_length = gene_info_df['exon_length']
    rc_per_gene = rc_df.iloc[:, -1]
    rc_total_protein_coding = rc_per_gene[protein_coding_gene_info_df.index].sum()
    rc_quantile_75th = np.quantile(rc_per_gene[protein_coding_gene_info_df.index], 0.75)

    with np.errstate(divide='ignore'):
        shared_part = np.log(rc_per_gene) + np.log(10**9) - np.log(gene_length)
        fpkm = np.exp(shared_part - np.log(rc_total_protein_coding))
        fpkm_uq = np.exp(shared_part - np.log(rc_quantile_75th))

    gene_basic_info_df = gene_info_df.loc[:, ['gene_name', 'hgnc_id']]
    fpkm_df = pd.concat(
        [gene_basic_info_df, rc_per_gene, fpkm, fpkm_uq], 
        axis=1, join='inner', sort=False
    )
    fpkm_df.columns = ['symbol', 'hgnc_id', 'read_count', 'fpkm', 'fpkm_uq']
    fpkm_df.index.name = 'gene_id'
    return fpkm_df.reset_index()


def add_symbol(df):
    mg = mygene.MyGeneInfo()
    genes = sorted(set([g.split('.')[0] for g in df['gene_id'].to_list()]))
    results = mg.querymany(genes, scopes='ensemblgene', fields='symbol', species='human')
    gene_id_to_symbol = {d['query']: d['symbol'] for d in results if 'symbol' in d}
    df['gene_symbol'] = [gene_id_to_symbol.get(g.split('.')[0], '')
                         for g in df['gene_id'].to_list()]
    df = df[['gene_id', 'gene_symbol', 'read_count', 'fpkm', 'fpkm_uq']]
    return df


def main(args):
    gene_info_df = pd.read_table(args.gene_info, index_col='gene_id')
    rc_df = pd.read_table(args.rc, comment='#', index_col='Geneid')
    fpkm_df = readcount_to_fpkm(gene_info_df, rc_df)

    # add gene symbol
    fpkm_df = add_symbol(fpkm_df)

    fpkm_df.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Readcount to FPKM and FPKM-UQ"
    )
    parser.add_argument('gene_info', help="Path to gene annotation info")
    parser.add_argument('rc', help="Path to readcount tsv")
    parser.add_argument('out', help="Path to output tsv")
    args = parser.parse_args()
    main(args)
