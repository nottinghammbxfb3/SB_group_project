# Franklin Brown Github repository
Generating structural and functional predictions for candidate genes generated by a selection scan on Cochlearia *excelsia* samples.

# Introduction
## The aims of this workflow:
- To create a consensus sequence of most common alleles at each nucleotide site from tetraploid and diploid samples for candidate genes.
- Generate mutated protein sequences.
- Identify position and substitution of mutations.
- Generate structure prediction models for mutated and reference proteins.
- Make functional predictions of mutations in different samples.

# Installation and tools

## To visualize and verify consensus sequence:
- IGV (version 2.17.2; Robinson et al., 2011)
  - https://igv.org/
 
## To create consensus sequence:
- Pandas (version 2.1.1; McKinney, W. & others, 2010, pp.51--56)
  - pip install pandas
  - https://pandas.pydata.org/
- Python (version 3.11.6)
  - https://www.python.org/
- Argparse (version 1.1)
  - pip install argparse
  - https://docs.python.org/3/library/argparse.html
 
## To create and view protein prediction models:
- AlphaFold2 (version 1.5.5)
  - (Mirdita et al., 2022)
- Phyre2 (version 2.0)
  - (Kelley et al., 2015)
- PyMOL (version 2.5.8)
  - https://pymol.org/

# Data
- Reference genome:
  - .fasta file containing all chromsosomes full sequence
- Reference annotation:
  - .ggf3 file containing full genome annotation
- Variant call files:
  - .vcf files for both tetraploid and diploid variants
- Gene names:
  - g30996 (chromosome 3)
  - g50328 (chromosome 6)
- Data source:
  - (Bray et al., 2023)
  - https://www.ebi.ac.uk/ena/browser/view/PRJEB66308
 
# Code
## To create consensus proteins:
- Using VCF_consensus.py:
```diff
  - chromosome_number  Chromosome number
  - gene_name          Gene name (Same format as in GFF file)
  - gff_path           Path to GFF file
  - vcf_path           Path to VCF file
  - reference_path     Path to reference FASTA file
  - snp_finder_output  Output file for SNP finder (.csv)
  - consensus_output   Output file for consensus sequence (.txt)
  - dna_to_aa_output   Output file for DNA to AA translation (.txt)
```
- Example usage:
```diff
python VCF_consensus.py 6 g50328 gff_file.gff3 vcf_file.vcf reference_genome.fasta snps_out.csv consensus_out.txt translation_out.txt
 ```
## VCF_consensus.py function explanation: (For detail on code used view VCF_consensus.py)
### SNP finder
  - Reads VCF file into pandas dataframe
  - Selects specified chromosome rows
  - Reads GFF file in to pandas dataframe
  - Selects specified chromsome, gene name and coding regions
  - Selects VCF rows corresponding to coding regions identified in GFF file
  - Create seperate VCF subsets of multiallelic and biallelic sites
  - Create alternative allele frequency column for both subsets
  - For biallelic subset, select rows above 0.5 frequency (most frequent)
  - For multiallelic subset, create reference allele frequency column:
    - Identify most frequent alternative allele
    - Remove rows where reference is more frequent than most frequent alternative allele
    - Checking throughout there are any alternative alleles remaining to avoid downstream error
  - Merge biallelic and multiallelic subsets and output as .csv for consensus function

  - Example output:

    | #CHROM | POS | REF | ALT |
    |--------|-----|-----|-----|
    | Cexcelsa_scaf_6 | 3691953 | T | G |
    | Cexcelsa_scaf_6 | 3692547 | A | C |
    | Cexcelsa_scaf_6 | 3693361 | A | C |
### Consensus
  - Read in reference fasta file for specified chromsome sequence as a string
  - Using the SNP finder output .csv file replace single nucleotide polymorphisms in reference chromosome
  - Using GFF file, as in SNP finder extract coding regions from chromosome
  - Output consensus sequence as a string, consisting of mutated coding regions, as .txt file
  - Example output:
```diff
    TCATGTAACAATGGAATCACCTGCTGAGCTGGTTTCGCCTTCGAGCCCAAAAA...
```
### DNA to amino acid
  - Read in Consensus output file
  - Translate first 3 reading frames using codon to amino acid translation dictionary
  - Create reverse complement of DNA sequence
  - Translate final 3 reading frames using codon to amino acid translation dictionary
  - Output to .txt file with each reading frame translation on newline
  - Example output:
```diff
  SCNNGITC_AGFAFFERAQKIGKFLRDGLICLC__I...AHRV
  HVTMESPAELVSPSSKEPKKLASF_ETV_FVFVNKS...TTAL
  M_QWNHLLSWFRLLRKSPKNWQVSERRFDLSLLINL...ARSS
+ MTSAVVSKPHRLKYDVFISFTGDTRHGFTEKLYKAL...TIV_
  _RAPSYRSRTGSNTTYSSASLETRATASRRSSIKHS...PLLH
  DERRRIEAAQAQIRRIHQLHWRHAPRLHGEAL_STR...HCYM
```
  - Infer correct reading frame following criteria (shown in green):
    - Start with start codon (M)
    - End with stop codoon (_)
    - Likely, fewest mismatches in multiple sequence alignment with reference
   
## Sequence alignment (For detail on code used view seq_align.py):
  - Run within text editor by copying sequences into funciton input
  - Algins 2 protein sequences of equal length
  - Outputs position, amount and mutations between 2 sequences to .txt file
  - Example output:
```diff
Number of mutations: 2

Pos: 9
Ref: M
Alt: T

Pos: 42
Ref: D
Alt: N
```
## Generating protein prediction models:
  - Using identified protein sequences for each sample, for each gene
  - Using web based AlphaFold2 colab and Phyre2 protein prediction software
  - Output: .pbd files, visualized using PyMOL

## PyMOL protein model visualization example code:
  - Load in .pbd files of models from same gene:
    - File > Open > protein_model.pbd
  - Align models:
    - align diploid_model.pbd, reference_model.pbd
  - Highlighting mutations identified with sequence alignment:
    - Colour mutated residues:
      - color red, (resi 7 or resi 21) and diploid_model.pbd
    - Colour truncated segments:
      - color blue, (resi 98-123) and reference_model.pbd







