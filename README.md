# Tile Miner: Checker Plot Generator for Nucleotide Sequences

**Author**: Jason Fajardo  
**License**: MIT  
**Version**: 1.1.0  
**Repository**:

---

## **Acknowledgments:**
This tool uses data and services provided by the National Center for Biotechnology Information (NCBI), including the Entrez and BLAST APIs. We acknowledge the invaluable role of NCBI in supporting open bioinformatics research.

For more information, see: https://www.ncbi.nlm.nih.gov

This project was developed in R and relies on the following packages:

- rentrez: for querying the NCBI Entrez API
- Biostrings: for DNA sequence manipulation and analysis
- ggplot2: for checker-style plotting
- reshape2

We acknowledge the authors and maintainers of these packages for their contributions to open-source science.

 **Note on tool development:**  
 The TileMiner program was conceptualized and directed by **Jason Fajardo**, with programming support provided through OpenAIs ChatGPT.  
 All R code was generated collaboratively using ChatGPT as a coding assistant, under Jason’s supervision, testing, and refinement.  
 This note is included for transparency and proper attribution of authorship roles.


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15693621.svg)](https://doi.org/10.5281/zenodo.15717685)



## Overview
Tile Miner is an R-based tool for generating checker-style tile plots from nucleotide sequences. Each sequence is fetched from NCBI using either an accession number or an Entrez query, and is visualized in a grid layout where each base is assigned a tile. Layout width and coloring affect the visual emergence of patterns, enabling exploration of sequence structure.

Version 1.1 introduces support for **coordinate slicing** (visualizing subregions of a sequence) and **user-defined color schemes**, allowing custom visual encoding strategies.

---

## Features
- Fetch nucleotide sequences by accession or search query.  
- Batch mode supports multiple sequence retrievals and layout generation.  
- Generates one `.png` plot per width per page.  
- Saves `.fasta` files for each sequence with gene names in the filename (if available).  
- Optional filters for sequence length.  
- Adjustable layout widths and pagination.  
- Supports sorting of search results.  
- **NEW**: Optional `start` and `end` parameters for slicing a specific region.  
- **NEW**: Customizable `color_scheme` argument for user-defined tile color encoding.  
- **NEW**: Per-sequence `.meta` files with metadata for downstream reuse.

---

## Usage

#### Single Accession Example

    tile_miner("NM_001101.5")


#### Batch Query Example

    tile_miner(
      query = "actin[Title] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
      widths = c(60, 80),
      max_results = 5,
      batch_size = 2,
      bases_per_page = 5000,
      min_length = 500,
      sort = "relevance"
    )


#### Other Examples

**Long RNA viral genome with multiple widths**
  
    tile_miner("NC_045512.2", widths = c(40, 60, 80), bases_per_page = 2000)

**Bacterial gene with genomic source and length limits**
  
    tile_miner("rpoB[Gene] AND Escherichia coli[Organism] AND biomol_genomic[PROP]", 
               widths = 70, max_results = 3, min_length = 5000)

**Intentional no-hit query to test error handling**

    tile_miner("unknown_gene[Gene] AND Homo sapiens[Organism]", max_results = 5)

**Slicing a genome region**

    tile_miner("NC_000913.3", start = 100000, end = 105000, widths = 50)

**Custom color scheme usage**

    tile_miner("NM_001101.5", widths = 60, color_scheme = c(A="#0000FF", T="#FF0000", G="#00FF00", C="#FFFF00"))

---

## Parameters
- `query`: Accession number or NCBI Entrez search string  
- `widths`: Integer vector; number of columns in each tile plot  
- `max_results`: Max number of results to retrieve (for search queries)  
- `batch_size`: Number of sequences to process per batch  
- `bases_per_page`: Number of bases per image page (controls pagination)  
- `save_fasta`: Save each sequence as a `.fasta` file (TRUE/FALSE)  
- `min_length`, `max_length`: Filter sequences by length (bp)  
- `verbose`: Print progress messages (TRUE/FALSE)  
- `sort`: Search result order (see below)  
- `start`, `end`: 1-based coordinates for region slicing (accession mode only)  
- `color_scheme`: Optional named vector mapping bases to hex color codes  

  -`sort` Options for `tile_fetch_batch()`

     - "relevance" (default)
     - "date"  
     - "submit_date"  
     - "sequence_length"  
     - "accn"  
     - "source"  

---

## Output
Each run creates a folder containing:  
- PNG images: One per width and page (e.g., NM_001101.5_w60_p1.png)  
- FASTA files: Named as Accession_Gene.fasta if gene metadata exists  
- `.meta` files: Per-sequence metadata stored in RDS format  
- `summary.txt` and `summary.csv`: Overview of plot settings and sources  

---

## NCBI Filter Grammar Cheatsheet
- `[Organism]`: Homo sapiens[Organism], Escherichia coli[Organism]  
- `[Gene]`: actin[Gene], p53[Gene]  
- `[Title]`: actin[Title] (less strict than `[Gene]`)  
- `[Sequence Length]`: 1000:3000[Sequence Length]  
- `[PDAT]`: 2020/01/01:2025/01/01[PDAT]  
- `[PROP]` Filters:  
  - `biomol_mrna[PROP]`: mRNA only  
  - `biomol_genomic[PROP]`: Genomic only  
  - `srcdb_refseq[PROP]`: RefSeq only  

You can combine filters using `AND`. Example:  

"p53[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP]"

---

## Visualization Note
Tile plots are highly dependent on:  
- The layout width (columns)  
- The tile colouring scheme  

Structural properties of sequences may become visible or obscured depending on how the sequence is folded into a grid. Custom color schemes enable encoding of biological signals such as GC content or entropy.

---

## License
MIT License. You are free to use, modify, and distribute this software. Attribution is appreciated.

---

## Citation
If you use this tool in a publication, please cite:  
Fajardo, J. (2025). *Tile Miner: A Visual Framework for Nucleotide Layout Analysis*. Zenodo. DOI: https://doi.org/10.5281/zenodo.15693621
