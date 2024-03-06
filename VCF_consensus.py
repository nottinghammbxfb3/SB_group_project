'''
Python script to create consensus protein sequence from VCF file.
Selects most prevelant allele at each position, avoiding IUPAC ambiguity codes.
Using reference fasta and gff annotation file.
Ouputs 6 open reading frame translation of consensus protein sequence.

Pandas and argparse need to be installed
    pip install pandas
    pip install argparse

For help:
   python VCF_consensus.py -h
'''
#!/usr/bin/env python3

import argparse
from argparse import RawTextHelpFormatter

def parse_arguments():
    parser = argparse.ArgumentParser(description="Takes VCF, GFF, reference genome and creates consensus protein from specified gene name and chromosome number, using most frequent bases in VCF file\nExample Usage:\npython VCF_consensus.py 6 g50328 gff_file.gff3 vcf_file.vcf reference_genome.fasta snps_out.csv consensus_out.txt translation_out.txt",formatter_class=RawTextHelpFormatter)
    parser.add_argument("chromosome_number", type=int, help="Chromosome number")
    parser.add_argument("gene_name", type=str, help="Gene name (Same format as in GFF file)")
    parser.add_argument("gff_path", type=str, help="Path to GFF file")
    parser.add_argument("vcf_path", type=str, help="Path to VCF file")
    parser.add_argument("reference_path", type=str, help="Path to reference FASTA file")
    parser.add_argument("snp_finder_output", type=str, help="Output file for SNP finder (.csv)")
    parser.add_argument("consensus_output", type=str, help="Output file for consensus sequence (.txt)")
    parser.add_argument("dna_to_aa_output", type=str, help="Output file for DNA to AA translation (.txt)")
    return parser.parse_args()


translation_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }


def snpfinder(gff_path,vcf_path,gene_name,chromosome_number,output_file):

    import pandas as pd

    # Function to save dataframe
    def save_df(df,filename):
        df.to_csv(filename, index = False)


    # Find the line number where '#CHROM' appears
    with open(vcf_path, "r") as vcf_file:
        for line_number, line in enumerate(vcf_file):
            if line.startswith("#CHROM"):
                break

    # Read in VCF, filter columns, chromosome and empty rows
    vcf = pd.read_csv(vcf_path, sep = '\t', skiprows = line_number)
    vcf = vcf[['#CHROM','POS','REF','ALT','QUAL','INFO']]
    vcf = vcf[vcf['#CHROM'] == 'Cexcelsa_scaf_' + str(chromosome_number)]
    vcf = vcf[vcf['ALT'] != '.']

    # Read in GFF, filter columns and chromsome
    gff = pd.read_csv(gff_path,sep='\t',header=None,
                  names=['CHROM','SOURCE','TYPE','START','END','SCORE','STRAND','PHASE','ID'])
    gff = gff[['CHROM','TYPE','START','END','ID']]
    gff = gff[gff['CHROM']=='Cexcelsa_scaf_'+str(chromosome_number)]

    # Filter GFF for gene and CDS
    gff = gff[gff['ID'].str.contains(gene_name)]
    gff = gff[gff['TYPE']=='CDS']

    # Select VCF rows from indetified GFF CDS regions
    cds_rows=[]
    for index, gff_row in gff.iterrows():
        cds_row = vcf[(vcf['POS'] >= gff_row['START']) & 
                      (vcf['POS'] <= gff_row['END'])]
        cds_rows.append(cds_row)
    vcf_cds = pd.concat(cds_rows)

    # Create Subset of multiallelic sites
    vcf_multi = vcf_cds[vcf_cds['ALT'].str.contains(',')]
    
    # Create subset of biallelic sites and remove empty rows
    vcf_bi = vcf_cds[~vcf_cds['ALT'].str.contains(',')]
    vcf_bi = vcf_bi[~vcf_bi['INFO'].str.contains('AC=0;')]

    # Make alt allele frequency column and select rows >0.5 frequency for bi subset
    vcf_bi_split = vcf_bi['INFO'].str.split(';')
    bi_af_vals=[]
    for row in vcf_bi_split:
        bi_af_val = None
        for item in row:
            if item.startswith('AF='):
                bi_af_val = float(item.split('=')[1])
                break
        bi_af_vals.append(bi_af_val)
    vcf_bi['AF'] = bi_af_vals
    vcf_bi = vcf_bi[vcf_bi['AF'] >= 0.5]
    bis = vcf_bi[['#CHROM','POS','REF','ALT']]

    # Make alt allele frequency column in multi subset
    vcf_multi_split = vcf_multi['INFO'].str.split(';')
    m_af_vals = []
    for row in vcf_multi_split:
        m_af_val = None
        for item in row:
            if item.startswith('AF='):
                m_af_val = item.split('=')[1]
                break
        m_af_vals.append(m_af_val)
    vcf_multi['AF'] = m_af_vals
    
    # Make sum of alt allele and reference frequency column and remove rows where reference is most frequent
    vcf_multi['AF'] = vcf_multi['AF'].str.split(',')
    vcf_multi['AF'] = vcf_multi['AF'].apply(lambda x: [float(value) for value in x])
    vcf_multi['AF_sum'] = vcf_multi['AF'].apply(lambda x: sum(x))
    vcf_multi['RF'] = 1 - vcf_multi['AF_sum']
    
    # Check if there are any remaining alt alleles
    if not vcf_multi['RF'].lt(0.5).any():
        if bis.empty:
            return 'No SNPs'
        else:
            save_df(bis,output_file)
    else:
        vcf_multi = vcf_multi[vcf_multi['RF'] <= 0.5]
    
    # Identify most frequent alt allele and select rows where alt is most frequent
    vcf_multi['AF_max'] = vcf_multi['AF'].apply(max)
    vcf_multi = vcf_multi[vcf_multi['AF_max'] >= vcf_multi['RF']]

    # Check if there are any remaining alt alleles
    if vcf_multi.empty:
        if bis.empty:
            return 'No SNPs'
        else:
            return save_df(bis,output_file)
    else:
        vcf_multi['max_i'] = vcf_multi['AF'].apply(lambda x: x.index(max(x)))
        vcf_multi['max_alt'] = vcf_multi.apply(lambda row: row['ALT'].split(',')[row['max_i']], axis=1)

        # Make dataframe containing all alt alleles
        multis = vcf_multi[['#CHROM','POS','REF','max_alt']]
        multis.rename(columns={'max_alt':'ALT'},inplace=True)
        
        all_alts = pd.concat([multis,bis])
        save_df(all_alts,output_file)


def consensus(snps_csv,chromosome_file,gff_file,chromosome_number,gene_name,output_file):
    import pandas as pd

    # Function to read fasta and extract given chromosome
    def get_chromosome_from_fasta(fasta_file, chromosome_number):
        header_name = 'Cexcelsa_scaf_'+str(chromosome_number)
        sequence = None
        with open(fasta_file, "r") as file:
            current_header = None
            for line in file:
                if line.startswith(">"):
                    current_header = line.strip()[1:]
                    if current_header == header_name:
                        sequence = ""
                elif current_header == header_name:
                    sequence += line.strip()
        return sequence
    chromosome = get_chromosome_from_fasta(chromosome_file, chromosome_number)

    # Function to replace snps at position
    def mutate(seq,position,ref,alt):
        mutated_seq = list(seq)
        if seq[position -1] == ref:
            mutated_seq[position -1] = alt
        return ''.join(mutated_seq)
        
    # Replace snps in chromosome
    snps = pd.read_csv(snps_csv)
    for index, row in snps.iterrows():
        position = row['POS']
        ref = row['REF']
        alt = row['ALT']
        chromosome = mutate(chromosome,position,ref,alt)

    # Read in GFF, filter columns and chromsome
    gff = pd.read_csv(gff_file,sep='\t',header=None,
                  names=['CHROM','SOURCE','TYPE','START','END','SCORE','STRAND','PHASE','ID'])
    gff = gff[['CHROM','TYPE','START','END','ID']]
    gff = gff[gff['CHROM']=='Cexcelsa_scaf_'+str(chromosome_number)]

    # Filter GFF for gene and CDS
    gff = gff[gff['ID'].str.contains(gene_name)]
    gff = gff[gff['TYPE']=='CDS']

    # Extract CDS regions from chromosome
    cds = []
    for index, row in gff.iterrows():
        start = row['START'] - 1
        end = row['END']
        cds.append(chromosome[start:end])
    cds = ''.join(cds)
    
    # Output consensus sequence to file
    with open(output_file, 'w')as out_f:
        out_f.write(cds)


def dna_to_aa(dna_file,translation_table,out_file):

    with open(dna_file, 'r')as file:
        dna = file.readline().strip()
    aas=[]

    # Translate first 3 reading frames
    for j in range(3):
        aa = ''
        for i in range(j, len(dna), 3):
            codon = dna[i:i + 3]
            if len(codon) == 3:
                aa += translation_table[codon]
        aas.append(aa)
    
    # Translate reverse compplement reading frames
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_comp = ''.join(complement[base] for base in reversed(dna))

    for j in range(3):
        aa = ''
        for i in range(j, len(rev_comp), 3):
            codon = rev_comp[i:i + 3]
            if len(codon) == 3:
                aa += translation_table[codon]
        aas.append(aa)

    with open(out_file, 'w') as f:
        for aa in aas:
            f.write(aa + '\n')


if __name__ == "__main__":
    args = parse_arguments()

    snpfinder(args.gff_path, args.vcf_path, args.gene_name, args.chromosome_number, args.snp_finder_output)
    consensus(args.snp_finder_output, args.reference_path, args.gff_path, args.chromosome_number, args.gene_name, args.consensus_output)
    dna_to_aa(args.consensus_output, translation_table, args.dna_to_aa_output)
