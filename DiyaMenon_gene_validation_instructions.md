# Gene Validation Pipeline for Biliary Atresia Candidates

## Instructions for Diya Menon

**Authors**: Ronnie Li and Stephen Hoang

**Created**: January 29, 2026

## Introduction

Welcome, Diya! This document provides detailed instructions for validating candidate genes that may contribute to biliary atresia. Steve and I identified several genes of interest, and your task is to systematically evaluate their biological plausibility as disease candidates.

Biliary atresia is a rare pediatric liver disease characterized by progressive destruction of bile ducts. Understanding the genetic contributors to this condition requires careful analysis of gene function, mutation patterns, sequencing quality, and genomic features. This pipeline will guide you through a comprehensive evaluation process.

## Note on Artificial Intelligence (AI)

Throughout this project, you are **encouraged to use AI assistants** (e.g., ChatGPT and Claude) to help you write code, understand concepts, and troubleshoot errors. However, it is **critical** that you:

- **Understand every line of code** you use, even if AI helped write it
- **Ask the AI to explain** any functions or syntax you don't recognize
- **Test code incrementally** rather than running large scripts all at once
- **Verify your results** make biological and computational sense

AI is a powerful tool for learning and productivity, but blindly copying code without understanding will prevent you from developing the skills you need. Think of AI as a teaching assistant. Use it to learn, not just to get answers.


## Project Overview

This pipeline consists of three main analytical components:

1. **Gene Annotation & Function** (Steps 1-2): Documenting basic gene properties and biological relevance
2. **Mutation Burden Analysis** (Steps 3-5): Calculating the ratio of observed mutations to constraint scores
3. **Sequencing Quality Metrics** (Steps 6-8): Assessing callability of gene regions
4. **Repetitive Element Analysis** (Steps 9-10): Quantifying repetitive sequences in genes

Each section includes a rationale explaining why we're performing these analyses.


## SECTION 1: Gene Annotation & Function

### Step 1: Record Basic Gene Information

**Task:** For each candidate gene, compile the following information:
- Gene name/symbol (e.g., CFTR)
- Ensembl ID (e.g., ENSG00000001626)
- Chromosome (e.g., chr7)
- Start position (genomic coordinate)
- Stop position (genomic coordinate)
- Gene biotype (protein_coding, lncRNA, pseudogene, etc.)

**Data Sources:**
- **Ensembl** (https://www.ensembl.org): Search by gene symbol
- **UCSC Genome Browser** (https://genome.ucsc.edu): Look up gene coordinates
- **NCBI Gene** (https://www.ncbi.nlm.nih.gov/gene): Alternative comprehensive resource

**Instructions:**
1. Create a spreadsheet or table with columns for each attribute listed above
2. For each gene, search Ensembl using the gene symbol
3. Record all requested information
4. Verify that coordinates correspond to the hg38/GRCh38 genome assembly

**Expected Output:** An Excel spreadsheet with one row per gene containing all basic annotation information.

---

### Step 2: Describe Gene Function and Liver Relevance

**Task:** Research and document the biological function of each gene, with specific emphasis on:
- The function of the protein product (or RNA product for non-coding genes)
- Relevance to liver biology, particularly in hepatocytes and cholangiocytes (bile duct cells)
- Any known associations with liver phenotypes or ciliary structures
- Evidence from literature, databases, and clinical resources

**Rationale:** Not all genes are equally plausible as biliary atresia candidates. Genes that are highly expressed in bile ducts, involved in cilia formation/function, or previously linked to liver diseases are stronger candidates. This step helps prioritize genes for further investigation.

**Data Sources:**
- **PubMed** (https://pubmed.ncbi.nlm.nih.gov): Search "[gene name] AND (liver OR cholangiocyte OR bile duct OR biliary OR cilia)"
- **ClinVar** (https://www.ncbi.nlm.nih.gov/clinvar): Check for pathogenic variants and associated phenotypes
- **OMIM** (https://www.omim.org): Look for gene-disease associations
- **GeneCards** (https://www.genecards.org): Comprehensive gene summaries
- **Human Protein Atlas** (https://www.proteinatlas.org): Check expression in liver tissues and cell types (hepatocytes, cholangiocytes).

**Instructions:**
1. For each gene, search PubMed for relevant publications
2. Check ClinVar for any reported pathogenic variants and their associated diseases
3. Review OMIM for known gene-disease relationships
4. Document findings in 2-3 paragraphs per gene
5. **Cite your sources** using PMID numbers, ClinVar accessions, or OMIM IDs

**Tips:**
- Focus on primary research articles, not just reviews
- Look for keywords: cholangiocyte, bile duct development, primary ciliary dyskinesia, cholestasis
- Note any animal model studies showing liver phenotypes

**Expected Output:** A written summary for each gene (1-2 paragraphs) with citations.


## SECTION 2: Mutation Burden Analysis

**Rationale:** Genes that are intolerant to loss-of-function (LOF) mutations in the general population (high constraint) but show many LOF mutations in biliary atresia patients may be disease-causing. The LOEUF (Loss-Of-function Observed/Expected Upper bound Fraction) score quantifies how tolerant a gene is to LOF mutations; lower scores indicate greater intolerance. By comparing observed mutations in patients to the gene's constraint score, we can identify genes where mutation burden is unexpectedly high.

### Step 3: Extract Mutation Counts from Formatted String

**Task:** Parse mutation annotation strings to extract counts of different mutation types.

**Background:** Your spreadsheet contains a column with mutation information formatted something like this:
```
(N1)Frameshift|(N2)LOF|(N3)Missense
```
For example: `(5)Frameshift|(3)LOF|(12)Missense`

Where:
- N1 = number of frameshift mutations
- N2 = number of loss-of-function mutations  
- N3 = number of missense mutations

**Instructions:**

You need to extract N1, N2, and N3 for each gene. Here's the general approach:

(*Note*: For all subsequent steps in this document, I will only provide R code examples, not Python code. But all tasks are possible in Python as well.)

**In R:**

```r
# Example approach - ask AI to explain each function
library(stringr)

# If your column is called "mutation_info"
# Extract the numbers using regular expressions
mutations_df <- your_data %>%
  mutate(
    frameshift = as.numeric(str_extract(mutation_info, "\\((\\d+)\\)Frameshift", group = 1)),
    lof = as.numeric(str_extract(mutation_info, "\\((\\d+)\\)LOF", group = 1)),
    missense = as.numeric(str_extract(mutation_info, "\\((\\d+)\\)Missense", group = 1))
  )
```

**In Python:**
```python
# Example approach - ask AI to explain each function
import re
import pandas as pd

# Extract numbers using regular expressions
def extract_mutations(mutation_string):
    frameshift = int(re.search(r'\((\d+)\)Frameshift', mutation_string).group(1))
    lof = int(re.search(r'\((\d+)\)LOF', mutation_string).group(1))
    missense = int(re.search(r'\((\d+)\)Missense', mutation_string).group(1))
    return frameshift, lof, missense

df[['frameshift', 'lof', 'missense']] = df['mutation_info'].apply(
    lambda x: pd.Series(extract_mutations(x))
)
```

**Learning Points:**
- **Regular expressions**: Patterns for finding text. `\d+` means "one or more digits"
- **Extraction**: Pulling specific parts of a string based on a pattern
- Ask AI: "Can you explain how this regular expression works step by step?"

**Expected Output:** Three new columns in your dataframe: frameshift, lof, and missense.

---

### Step 4: Obtain LOEUF Scores

**Task:** Download constraint metrics from gnomAD and extract LOEUF scores for your genes.

**Background:** The LOEUF score comes from the Genome Aggregation Database (gnomAD), which analyzed genetic variation in hundreds of thousands of individuals. Genes with low LOEUF scores (<0.35) are highly constrained and intolerant to LOF mutations.

**Instructions:**

1. **Download the data:**
   - Go to UCSC Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables
   - Set the following options:
     - **clade:** Mammal
     - **genome:** Human
     - **assembly:** GRCh38/hg38
     - **group:** Variation
     - **track:** gnomAD Constraint Metrics
     - **table:** Gene LoF (pliByGene)
     - **output format:** all fields from selected table
   - Click "get output"
   - Save the file (e.g., `gnomad_loeuf.txt`)

2. **Load into R:**
```r
# Read the downloaded file
loeuf_data <- read.table("gnomad_loeuf.txt", header = TRUE, sep = "\t")

# Merge with your gene list
# Make sure gene symbols match between datasets
merged_data <- left_join(your_genes, loeuf_data, by = c("gene_symbol" = "gene_name"))

# Handle missing values
# Option 1: Keep as NA
merged_data$oeuf_upper <- ifelse(is.na(merged_data$oeuf_upper), NA, merged_data$oeuf_upper)

# Option 2: Replace with median
median_loeuf <- median(loeuf_data$oeuf_upper, na.rm = TRUE)
merged_data$oeuf_upper <- ifelse(is.na(merged_data$oeuf_upper), median_loeuf, merged_data$oeuf_upper)
```

**Important Notes:**
- The LOEUF column may be named `oeuf_upper` or `oe_lof_upper` depending on the source
- Gene symbols must match exactly (case-sensitive)
- Some genes may not have LOEUF scores (very rare genes or genes with no coding sequence)

**Expected Output:** Your dataframe now includes a LOEUF score column for each gene.

---

### Step 5: Calculate Mutation Burden Ratio

**Task:** Calculate the ratio of observed LOF mutations to the LOEUF score.

**Rationale:** This ratio tells us whether we're seeing more LOF mutations than expected given the gene's tolerance. A high ratio suggests the gene may be involved in disease.

**Instructions:**

Calculate two metrics:

1. **LOF/LOEUF ratio:** Number of observed LOF mutations divided by LOEUF score
2. **Total mutations/LOEUF ratio:** Total mutations (frameshift + LOF + missense) divided by LOEUF score

**In R:**
```r
mutation_burden <- merged_data %>%
  mutate(
    lof_to_loeuf = lof / oeuf_upper,
    total_to_loeuf = (frameshift + lof + missense) / oeuf_upper
  )

# Genes with high ratios are candidates of interest
high_burden_genes <- mutation_burden %>%
  filter(lof_to_loeuf > threshold) %>%  # You'll need to determine an appropriate threshold
  arrange(desc(lof_to_loeuf))
```

**Interpretation:**
- Higher ratios indicate more observed mutations relative to constraint
- Genes with LOEUF < 0.35 and high mutation counts are particularly interesting
- Consider both absolute mutation counts and ratios

**Expected Output:** Two new calculated columns showing mutation burden ratios.

---

## SECTION 3: Sequencing Quality Assessment (Callability)

**Rationale:** Not all parts of the genome are equally easy to sequence. Some regions have low coverage, making it difficult to confidently call variants. If a gene has poor "callability" (many bases with insufficient sequencing depth), we might be missing mutations simply because we couldn't detect them. Quantifying callability helps us understand whether low mutation counts reflect true biology or technical limitations.

### Step 6: Import Coding Sequence Data

**Task:** Load and examine the coding sequence (CDS) information for each gene.

**Background:** We will provide you with:
1. An Excel file with the total CDS length for each gene (canonical transcript only)
2. A file with genomic coordinates for each coding exon

**Instructions:**

**In R:**
```r
# Install required package if needed
# install.packages("readxl")
library(readxl)

# Load CDS lengths
cds_lengths <- read_excel("gene_cds_lengths.xlsx")

# Load CDS coordinates (format: gene, chr, start, end)
cds_coords <- read.table("gene_cds_coordinates.bed", header = TRUE, sep = "\t")
# or if it's a BED file without header:
cds_coords <- read.table("gene_cds_coordinates.bed", header = FALSE, sep = "\t",
                         col.names = c("chr", "start", "end", "gene"))

# Examine the data
head(cds_lengths)
head(cds_coords)

# Check: do genes have multiple CDS regions?
cds_coords %>% 
  gr by(gene_name) %>%
  summarize(num_exons = n())
```

**What to look for:**
- Each gene should have one entry in the CDS lengths file
- Each gene will have multiple entries in the coordinates file (one per coding exon)
- Coordinates should be in hg38 assembly

**Expected Output:** Loaded dataframes; understanding of data structure.

---

### Step 7: Calculate Callable Bases per Gene

**Task:** Determine how many bases in each gene's coding sequence meet quality thresholds.

**Background:** We will provide a BED file of "callable" regions—genomic intervals where at least 90% of samples (1,575 out of 1,750) have sequencing depth ≥10X. These are high-confidence regions where we trust variant calls.

**Instructions:**

You need to find the overlap between:
- CDS coordinates for each gene (from Step 6)
- Callable regions (provided BED file)

**Conceptual approach:**
1. For each gene, get all its CDS intervals
2. Find which parts of those intervals overlap with callable regions
3. Sum the total callable bases

**In R using GenomicRanges:**
```r
# Install if needed
# BiocManager::install("GenomicRanges")
library(GenomicRanges)

# Convert CDS coordinates to GRanges object
cds_gr <- GRanges(
  seqnames = cds_coords$chr,
  ranges = IRanges(start = cds_coords$start, end = cds_coords$end),
  gene = cds_coords$gene
)

# Load callable regions
callable_bed <- read.table("callable_regions.bed", sep = "\t",
                           col.names = c("chr", "start", "end"))
callable_gr <- GRanges(
  seqnames = callable_bed$chr,
  ranges = IRanges(start = callable_bed$start, end = callable_bed$end)
)

# Find overlaps
overlaps <- findOverlaps(cds_gr, callable_gr)

# Calculate callable bases per gene
# This is the tricky part - ask AI for help with the details
callable_by_gene <- cds_gr[queryHits(overlaps)] %>%
  intersect(callable_gr[subjectHits(overlaps)]) %>%
  # Group by gene and sum widths
  # ... (you will need to complete this)
```

**Note:** If you're unsure about the graphical representation of interval overlaps, ask Ronnie for clarification.

**Expected Output:** A table with each gene and its total number of callable bases.

---

### Step 8: Calculate Percent Callable

**Task:** Determine what percentage of each gene's coding sequence is callable.

**Instructions:**

**In R:**
```r
# Merge callable bases with CDS lengths
callability <- left_join(cds_lengths, callable_df, by = "gene")

# Calculate percentage
callability <- callability %>%
  mutate(
    percent_callable = (callable_bases / cds_length) * 100
  )

# Identify genes with poor callability
low_callability <- callability %>%
  filter(percent_callable < 80) %>%
  arrange(percent_callable)
```

**Interpretation:**
- Genes with <80% callability may have unreliable mutation counts
- Low callability can explain why some genes appear to have few mutations
- Consider flagging genes with poor callability in your final results

**Expected Output:** Percent callable for each gene; list of genes with callability concerns.

---

## SECTION 4: Repetitive Element Analysis

**Rationale:** Repetitive sequences (transposons, tandem repeats, low-complexity regions) are difficult to sequence and align accurately. Genes with high repetitive content may have unreliable variant calls, similar to the callability issue. Additionally, some repetitive elements are genuinely functional and could be relevant to disease. Quantifying repetitive content helps us interpret mutation data and understand gene structure.

### Step 9: Download RepeatMasker Annotations

**Task:** Obtain a comprehensive list of all repetitive elements in the human genome.

**Instructions:**

1. **Download from UCSC Table Browser:**
   - Go to: https://genome.ucsc.edu/cgi-bin/hgTables
   - Set the following options:
     - **clade:** Mammal
     - **genome:** Human
     - **assembly:** GRCh38/hg38
     - **group:** Repeats
     - **track:** RepeatMasker
     - **table:** rmsk
     - **region:** genome (to get all repeats)
     - **output format:** BED
   - Click "get output"
   - Save as `repeatmasker_hg38.bed`

2. **Load into R:**
```r
# Load RepeatMasker BED file
# BED format: chr, start, end, name, score, strand, ...
repeats <- read.table("repeatmasker_hg38.bed", sep = "\t", header = FALSE)

# Assign column names if not provided
colnames(repeats) <- c("chr", "start", "end", "repeat_name", "score", "strand")

# Examine the data
head(repeats)
summary(repeats)

# How many different types of repeats?
table(repeats$repeat_name)
```

**What you're looking at:**
- This file contains millions of intervals representing repetitive DNA
- Common types: LINE, SINE, LTR (transposable elements), Simple_repeat, Low_complexity
- Each interval shows where in the genome a repeat is located

**Expected Output:** Loaded dataframe with repetitive element coordinates.

---

### Step 10: Calculate Percent Repetitive per Gene

**Task:** Determine what percentage of each gene's coding sequence overlaps with repetitive elements.

**Instructions:**

Similar to Step 7, you need to find overlaps between gene CDS regions and RepeatMasker regions.

**In R using GenomicRanges:**
```r
# Convert repeats to GRanges
repeats_gr <- GRanges(
  seqnames = repeats$chr,
  ranges = IRanges(start = repeats$start, end = repeats$end)
)

# Find overlaps with CDS (cds_gr from Step 7)
repeat_overlaps <- findOverlaps(cds_gr, repeats_gr)

# Calculate repetitive bases per gene
# Similar logic to Step 7
repetitive_by_gene <- cds_gr[queryHits(repeat_overlaps)] %>%
  intersect(repeats_gr[subjectHits(repeat_overlaps)]) %>%
  # Sum by gene
  # ... (you will need to complete)

# Calculate percentage
gene_repetitiveness <- cds_lengths %>%
  left_join(repetitive_by_gene, by = "gene") %>%
  mutate(
    percent_repetitive = (repetitive_bases / cds_length) * 100,
    percent_repetitive = ifelse(is.na(percent_repetitive), 0, percent_repetitive)
  )
```

**Interpretation:**
- Genes with >30% repetitive content may have reduced sequencing reliability
- Some genes are naturally repeat-rich (e.g., mucins, certain structural proteins)
- High repetitive content doesn't disqualify a gene but should be noted

**Expected Output:** Percent repetitive for each gene.

---

## Final Integration and Summary

After completing all steps, you should combine your results into a comprehensive summary table:

**Final table should include:**
- Gene symbol, Ensembl ID, chromosome, coordinates, biotype (Step 1)
- Mutation counts: frameshift, LOF, missense (Step 3)
- LOEUF score (Step 4)
- LOF/LOEUF ratio, Total/LOEUF ratio (Step 5)
- CDS length, callable bases, percent callable (Steps 6-8)
- Repetitive bases, percent repetitive (Steps 9-10)
- Written functional summary with citations (Step 2)

**In R:**
```r
# Combine all analyses
final_results <- cds_lengths %>%
  left_join(mutation_burden, by = "gene") %>%
  left_join(callability, by = "gene") %>%
  left_join(gene_repetitiveness, by = "gene")

# Export to Excel for review
library(writexl)
write_xlsx(final_results, "biliary_atresia_gene_validation_results.xlsx")
```
---

## Tips for Success

1. **Work incrementally:** Complete one step at a time and verify your output before moving on
2. **Check your data:** Always use `head()`, `summary()`, `str()` to inspect dataframes
3. **Handle errors gracefully:** When code fails, read the error message carefully
4. **Ask questions:** If something doesn't make sense, ask Ronnie or Steve, or use AI to clarify
5. **Document your work:** Add comments to your code explaining what each section does
6. **Keep a log:** Note any decisions you make (e.g., how you handled missing LOEUF scores)

## Common Issues and Solutions

**Problem:** Gene symbols don't match between datasets
- **Solution:** Check for case sensitivity, use `toupper()` in R or `.str.upper()` in Python to standardize

**Problem:** BED files have different coordinate systems (0-based vs 1-based)
- **Solution:** UCSC BED files are 0-based, half-open. If using with 1-based tools, be careful

**Problem:** Missing data in joins
- **Solution:** Use `left_join()` in R or `merge(..., how='left')` in Python to keep all your genes

**Problem:** Genomic intervals aren't merging correctly
- **Solution:** Ensure chromosome names match exactly (chr1 vs 1), sort by coordinates

**Problem:** AI generates code that doesn't work
- **Solution:** Ask AI to explain the error, break code into smaller pieces, test each part manually yourself.

---

## Additional Resources

### R Resources
- GenomicRanges tutorial: https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
- dplyr cheat sheet: https://rstudio.github.io/cheatsheets/data-transformation.pdf

### Python Resources
- pandas documentation: https://pandas.pydata.org/docs/
- pybedtools tutorial: https://daler.github.io/pybedtools/

### Genomics Resources
- UCSC Genome Browser tutorial: https://genome.ucsc.edu/training/
- BED format specification: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

---

## Questions to Ask Yourself

As you complete this pipeline, reflect on these questions:

1. Which genes have the strongest evidence for involvement in biliary atresia? Why?
2. Are there genes with high mutation burden but poor callability? How should we interpret those?
3. Do any genes show unexpectedly high repetitive content? What might that mean?
4. Which genes would you prioritize for follow-up experiments? What criteria are you using?
5. What are the limitations of this analysis? What additional data would be helpful?

---

## Conclusion

This pipeline will give you hands-on experience with real bioinformatics workflows: data integration, genomic interval operations, quality control, and biological interpretation. Take your time, understand each step, and don't hesitate to ask for help. Good luck!

**Remember:** The goal isn't just to complete the analysis; it's to understand the biology and develop the computational skills you'll use throughout your budding bioinformatics career.
