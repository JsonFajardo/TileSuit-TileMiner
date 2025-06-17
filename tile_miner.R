# === tile_miner(): Simple Tile Plot Generator for Nucleotide Sequences ===
# NOTE: The coloring scheme and column layout (width) significantly affect visual patterns.
#       Some structural properties may emerge or become hidden depending on tile layout.
# Description:
# This tool fetches nucleotide sequences from NCBI using either an accession number or a query string,
# and generates checker-style tile plots of the sequences. Each base is rendered as a tile in a grid
# where rows and columns represent sequence position.
# 
# There are two main modes:
# - `tile_fetch_one()`: For a single accession (e.g. NM_001101.5)
# - `tile_fetch_batch()`: For Entrez queries (e.g. "actin[gene] AND Homo sapiens[Organism]")
# 
# The wrapper function `tile_miner()` auto-selects between them based on the input.
# 
# Features:
# - Visual layout variations via the `widths` parameter
# - FASTA export (optional)
# - Minimum and maximum length filters
# - Pagination of long sequences
# - Metadata extraction (e.g., gene name, organism) from FASTA headers
# 
# Files generated:
# - One PNG image per sequence per layout width per page
# - Optional: one FASTA per sequence
# - `summary.txt` and `summary.csv` for each run
# 
# Example (Single Accession):
#   tile_miner("NM_001101.5")
# 
# Example (Batch Search):
#   tile_miner("actin[Title] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
#              widths = c(60), max_results = 5, min_length = 500)
#
# === Extended Test Examples ===
# These test a wide range of parameters and scenarios:
#
# 1. Single accession, large sequence, multiple widths:
#   tile_miner("NC_045512.2", widths = c(40, 60, 80), bases_per_page = 2000)
#
# 2. Single accession, with FASTA disabled and custom length limits:
#   tile_miner("NM_001101.5", widths = 50, save_fasta = FALSE, min_length = 1000, max_length = 3000)
#
# 3. Batch query, max 3 results, small batch size, two widths:
#   tile_miner("insulin[Title] AND biomol_mrna[PROP]", max_results = 3, batch_size = 1, widths = c(40, 80))
#
# 4. Longer query with organism, source, and publication date filter:
#   tile_miner("p53[Gene] AND Homo sapiens[Organism] AND srcdb_refseq[PROP] AND 2010/01/01:2024/12/31[PDAT]",
#              widths = c(60), max_results = 5, min_length = 1000)
#
# 5. Intentionally too-strict query (tests graceful failure):
#   tile_miner("unknown_gene_name[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
#              max_results = 5)
#
# 6. Bacterial gene, medium size, with genomic source filter:
#   tile_miner("rpoB[Gene] AND Escherichia coli[Organism] AND biomol_genomic[PROP]",
#              widths = 70, max_results = 3, min_length = 5000)
#
# 7. RNA viral sequence, small genome:
#   tile_miner("NC_001802.1", widths = c(40, 60, 100), bases_per_page = 1500)

# === NCBI Filter Grammar Cheatsheet ===
# These are common filters used with Entrez queries for `tile_fetch_batch()`:
#   - [Organism]                 => Homo sapiens[Organism], Escherichia coli[Organism]
#   - [Gene]                     => actin[Gene], p53[Gene]
#   - [Title]                    => actin[Title] (less strict than [Gene])
#   - [Sequence Length]          => 100:2000[Sequence Length]
#   - [PDAT]                     => 2020/01/01:2024/12/31[PDAT]
#   - [PROP] filters:
#       biomol_mrna[PROP]       => Only mRNA sequences
#       biomol_genomic[PROP]    => Genomic sequences
#       biomol_rrna[PROP]       => rRNA sequences
#       srcdb_refseq[PROP]      => From RefSeq
# 
# Example combinations:
#   - "actin[Title] AND Homo sapiens[Organism] AND biomol_mrna[PROP]"
#   - "p53[Gene] AND Mus musculus[Organism] AND srcdb_refseq[PROP]"
#   - "insulin[Title] AND 1000:3000[Sequence Length]"
# 
# Note: Using multiple strict filters (e.g., gene+organism+refseq) may lead to no hits. Relax filters gradually.
# 
# Required packages
library(rentrez)
library(Biostrings)
library(ggplot2)
library(reshape2)

# === Helper Functions ===
is_accession <- function(x) {
  grepl("^[A-Z]{1,2}_[0-9]+(\\.[0-9]+)?$", x) && !grepl(" ", x)
}

clean_seq <- function(seq) {
  toupper(gsub("[^ATGCUN]", "", seq))
}

save_fasta_file <- function(id, gene, full_name, seq, outdir) {
  label <- if (!is.na(gene) && nzchar(gene)) paste0(id, "_", gene) else id
  fasta_text <- paste0(">", label, " ", full_name, "
", seq)
  writeLines(fasta_text, file.path(outdir, paste0(label, ".fasta")))
}

to_matrix_df <- function(chars, width) {
  len <- length(chars)
  height <- ceiling(len / width)
  chars <- c(chars, rep(NA, height * width - len))
  mat <- matrix(chars, nrow = height, byrow = TRUE)
  df <- melt(mat, varnames = c("row", "col"), value.name = "char")
  df$row_offset <- max(df$row) - df$row + 1
  df$col_offset <- df$col
  df <- df[!is.na(df$char), ]
  df
}

render_tile_plot <- function(df, width, page, accession, outdir, seq_len, total_pages, full_name = NULL) {
  title <- paste0("Accession: ", accession, if (!is.null(full_name)) paste0(" | ", full_name) else "")
  subtitle <- paste("Width:", width, "| Page:", page, "of", total_pages)
  caption <- paste0("Length: ", seq_len, " bp
Total Pages: ", total_pages, "
Columns: ", width)
  g <- ggplot(df, aes(x = col_offset, y = row_offset)) +
    geom_tile(aes(fill = char), colour = "white", linewidth = 0) +
    scale_fill_manual(values = c(A="black", T="gray30", C="gray30", G="black", N="orange"), na.value = "white") +
    theme_void() +
    coord_fixed() +
    ggtitle(title, subtitle = subtitle) +
    labs(caption = caption) +
    theme(plot.caption = element_text(hjust = 0, size = 8),
          plot.subtitle = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold"))
  ggsave(file.path(outdir, paste0(gsub("[[:punct:]]", "_", accession),
                                  "_w", width, "_p", page, ".png")),
         g, width = 10, height = 10)
}

process_accession <- function(id, outdir, widths, bases_per_page, min_length, max_length, save_fasta, verbose) {
  metadata_summary <- list()
  if (verbose) message("Fetching accession: ", id)
  fasta <- entrez_fetch(db = "nuccore", id = id, rettype = "fasta", retmode = "text")
  lines <- strsplit(fasta, "\n")[[1]]
  header <- gsub("^>", "", lines[1])
  parts <- strsplit(header, " ", fixed = TRUE)[[1]]
  accession <- sub(" .*", "", header)
  full_name <- paste(parts[-1], collapse = " ")
  organism <- ifelse(grepl("\\[.*\\]", header), sub(".*\\[(.*)\\].*", "\\1", header), NA)
  gene <- ifelse(grepl("gene", header, ignore.case = TRUE), sub(".*gene[:= ]*([^\\s]+).*", "\\1", header, ignore.case = TRUE), NA)
  header_metadata <- paste0(full_name, if (!is.na(gene)) paste0(" | Gene: ", gene) else "", if (!is.na(organism)) paste0(" | Organism: ", organism) else "")
  seq <- clean_seq(paste(lines[-1], collapse = ""))
  seq_len <- nchar(seq)
  if (!is.null(min_length) && seq_len < min_length) {
    if (verbose) message("Skipping ", id, ": sequence length ", seq_len, " < min_length ", min_length)
    return(metadata_summary)
  }
  if (!is.null(max_length) && seq_len > max_length) {
    if (verbose) message("Skipping ", id, ": sequence length ", seq_len, " > max_length ", max_length)
    return(metadata_summary)
  }
  if (seq_len > 500000 && verbose) {
    warning("Sequence ", id, " is very long (", seq_len, " bp). Plotting may be slow or memory-intensive.")
  }
  if (save_fasta) save_fasta_file(id, gene, full_name, seq, outdir)
  chars <- strsplit(seq, "")[[1]]
  for (w in widths) {
    total_len <- length(chars)
    pages <- ceiling(total_len / bases_per_page)
    entry <- list(
      accession = accession,
      name = header_metadata,
      length = total_len,
      width = w,
      pages = pages
    )
    metadata_summary[[length(metadata_summary) + 1]] <- entry
    for (p in seq_len(pages)) {
      if (verbose) message("Plotting ", id, " | Width ", w, " | Page ", p)
      start <- ((p - 1) * bases_per_page) + 1
      end <- min(start + bases_per_page - 1, total_len)
      df <- to_matrix_df(chars[start:end], w)
      render_tile_plot(df, w, p, accession, outdir, total_len, pages, full_name = header_metadata)
    }
  }
  return(metadata_summary)
}

# === tile_fetch_one(): Process a single accession ===
tile_fetch_one <- function(accession,
                           widths = c(60), # number of columns in plot
                           bases_per_page = 5000,
                           min_length = 200,
                           max_length = NULL,
                           save_fasta = TRUE,
                           verbose = TRUE) {
  if (verbose) message("Detected accession input — skipping search.")
  outdir <- paste0(gsub("[[:punct:]]", "_", accession), "_tileminer")
  dir.create(outdir, showWarnings = FALSE)
  metadata_summary <- process_accession(accession, outdir, widths, bases_per_page, min_length, max_length, save_fasta, verbose)
  summary_file <- file.path(outdir, "summary.txt")
  summary_file_csv <- file.path(outdir, "summary.csv")
  cat("Tile Miner Output Summary\n==========================\n", file = summary_file)
  csv_rows <- list()
  for (i in seq_along(metadata_summary)) {
    data <- metadata_summary[[i]]
    cat(paste0(
      "Accession: ", data$accession, "\n",
      "Name: ", data$name, "\n",
      "Length: ", data$length, " bp\n",
      "Width: ", data$width, "\n",
      "Pages: ", data$pages, "\n\n"
    ), file = summary_file, append = TRUE)
    csv_rows[[i]] <- data
  }
  df_summary <- do.call(rbind, lapply(csv_rows, as.data.frame))
  write.csv(df_summary, summary_file_csv, row.names = FALSE)
}

# === tile_fetch_batch(): Perform search and batch fetch ===
# - `sort`: Optional parameter for controlling the order of search results.
#           Default is "relevance". Other accepted values:
#           "date", "submit_date", "sequence_length", "accn", "source"
#           Example: sort = "sequence_length" to prioritize longest sequences.
# - `max_results`: Maximum number of sequence IDs to retrieve from NCBI.
#                  This limits how many total entries will be downloaded (e.g., 100).
# - `batch_size`: Number of sequences to process at once.
#                 Useful to avoid API timeouts, rate limits, or memory issues.
#                 For example, setting `batch_size = 10` with `max_results = 100`
#                 processes 10 batches of 10 sequences each.
tile_fetch_batch <- function(query,
                             sort = "relevance",
                             widths = c(60),
                             max_results = 10,
                             batch_size = 5,
                             bases_per_page = 5000,
                             min_length = 200,
                             max_length = NULL,
                             save_fasta = TRUE,
                             verbose = TRUE) {
  if (verbose) message("Performing NCBI search for: ", query)
  search <- entrez_search(db = "nuccore", term = query, retmax = max_results, sort= sort)
  ids <- search$ids
  if (length(ids) == 0) {
    stop("No sequences found for query.")
  }
  query_name <- gsub("[[:punct:] ]", "_", substr(query, 1, 40))
  outdir <- paste0(query_name, "_tileminer")
  dir.create(outdir, showWarnings = FALSE)
  metadata_summary <- list()  # <-- Add this line here
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  for (batch in batches) {
    for (id in batch) {
      try({
        result <- process_accession(id, outdir, widths, bases_per_page, min_length, max_length, save_fasta, verbose)
        if (length(result) > 0) {
          metadata_summary <- c(metadata_summary, result)
        }
      }, silent = TRUE)
    }
  }
  # Write summary after batch processing
  summary_file_txt <- file.path(outdir, "summary.txt")
  summary_file_csv <- file.path(outdir, "summary.csv")
  cat("Tile Miner Output Summary
==========================
", file = summary_file_txt)
  csv_rows <- list()
  for (i in seq_along(metadata_summary)) {
    data <- metadata_summary[[i]]
    cat(paste0(
      "Accession: ", data$accession, "
",
"Name: ", data$name, "
",
"Length: ", data$length, " bp
",
"Width: ", data$width, "
",
"Pages: ", data$pages, "

"
    ), file = summary_file_txt, append = TRUE)
    csv_rows[[i]] <- data
  }
  df_summary <- do.call(rbind, lapply(csv_rows, as.data.frame))
  write.csv(df_summary, summary_file_csv, row.names = FALSE)
}

# === Wrapper Function ===
# - Supports all parameters from tile_fetch_one() and tile_fetch_batch(), including:
#     widths, max_results, batch_size, bases_per_page, min_length, max_length,
#     save_fasta, verbose, and sort (controls NCBI search order).
#     Valid sort values: "relevance" (default), "date", "submit_date", "sequence_length", "accn", "source"
tile_miner <- function(query, ...) {
  if (is_accession(query)) {
    tile_fetch_one(query, ...)
  } else {
    tile_fetch_batch(query, ...)
  }
}




# === Extended Test Examples ===
# These test a wide range of parameters and scenarios:

# 1. Single accession, large sequence, multiple widths:
   tile_miner("NC_045512.2", widths = c(40, 60, 80), bases_per_page = 2000)

# 2. Single accession, with FASTA disabled and custom length limits:
   tile_miner("NM_001101.5", widths = 50, save_fasta = FALSE, min_length = 1000, max_length = 3000)

# 3. Batch query, max 3 results, small batch size, two widths:
   tile_miner("insulin[Title] AND biomol_mrna[PROP]", max_results = 3, batch_size = 1, widths = c(40, 80))

# 4. Longer query with organism, source, and publication date filter:
   tile_miner("p53[Gene] AND Homo sapiens[Organism] AND srcdb_refseq[PROP] AND 2010/01/01:2024/12/31[PDAT]",
              widths = c(60), max_results = 5, min_length = 1000)

# 5. Intentionally too-strict query (tests graceful failure):
   tile_miner("unknown_gene_name[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
              max_results = 5)

# 6. Bacterial gene, medium size, with genomic source filter:
   tile_miner("rpoB[Gene] AND Escherichia coli[Organism] AND biomol_genomic[PROP]",
              widths = 70, max_results = 30, min_length = 5000)

 # 7. RNA viral sequence, small genome:
  tile_miner("NC_001802.1")

 # 8. Mitochondrial sequence.
  tile_fetch_batch("mitochondrion[Filter] AND Homo sapiens[Organism]")
 
  # 9. Multiple Sequences One Call, Three Plots. 
  
 tile_miner(
   query = "NM_001101.5 OR NM_007393 OR NM_131031.2",
   widths = c(60),
   max_results = 3,
   batch_size = 3
 )
 
 # 10. Batch Search.
 
 tile_miner(
   query = "actin beta[Title] AND (Homo sapiens[Organism] OR Mus musculus[Organism] OR Danio rerio[Organism]) AND biomol_mrna[PROP] AND refseq[filter]",
   widths = c(60),
   max_results = 6,
   batch_size = 2,
   min_length = 500
 )
 