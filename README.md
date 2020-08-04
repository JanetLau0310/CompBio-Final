# CompBio-Final
Gene expression analysis, PPI networks, and COVID-19 sequence data Codon Analysis.

## Project One: Classifying expression data.
For this problem we will work with a public gene expression data set. To create this data set, RNA was extracted from bronchoscopy samples, and Affymetrix microarrays were used to measure gene expression. The resulting gene expression levels were normalized, and the data set was deposited in the NCBI GEO repository with accession number GSE994. <br> For this assignment, we collapsed the measurements to gene symbols (i.e. if a gene was represented by multiple Affymetrix probe sets on the array, we assigned the median expression value to the gene). Further, we are looking at data from only a subset of the patients to make the assignment more straightforward. The expression values are reported in log-transformed "average difference" values.












## Project Two: Protein-protein interaction networks.
For this problem, you will be performing interactive data analysis rather than implementing an algorithm. You will use a data file representing the yeast protein-protein interaction (PPI) network, posted as yeast.ppi on the classmaterials web site.<br>
This tab-delimited text file includes 76,025 yeast physical protein-protein interactions among 5,001 proteins. Each row corresponds to one interaction and has three columns:
the first and the second columns are the two interacting proteins, while the third entry is an number in (0, 1] that shows the confidence of the existence of the interaction











## Project Three: COVID-19 sequence analysis.
For this problem, you will work with sequences from COVID-19 virus samples collected from January through March of this year. The data were obtained from the NCBI
SARS-Cov-2 nucleotide protein sequence site linked from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/ on March 31st, 2020.<br>
Data files including:
- covid19sequences.fasta (whole genome sequences for 320 genomes)
- covid19codingClean.fasta (gene coding nucleotide seqs for genes from the above
sequences)
- metadata.table.csv (metadata about where and when the sequences were col-
lected)
- ncov-sequences-info.txt (additional metadata including geographic state for US
and China)
