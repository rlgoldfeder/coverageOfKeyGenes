# Coverage of Key Genes #

This pipeline calculates and plots the coverage of key genes for your input sample. 

### Expected Output: ###

This pipeline outputs a table and pdf with information about the number (and %) of bases per gene failing to meet the coverage minimums.

### User inputs: ##

**Relevant coverage thresholds**
- Minimum base quality
- Minimum mapping quality
- Minimum depth of coverage

**Relevant Genes**

Relevant genes must be input a BED file, where the first 3 columns are chrom, start, and stop and the fourth column is the gene name.
The BED file can contain several lines for each genes (ie: one line per exon). 

**Sequence Data**

You must include a BAM file and its index (BAI) file.

**GATK jar File**
You may download GATK from the following link:

https://www.broadinstitute.org/gatk/download/

**Reference Genome**

A reference genome (used by GATK). This must match with the reference genome you used to create your BAM file.

**GATK memory requirements**

Default is 16G.
