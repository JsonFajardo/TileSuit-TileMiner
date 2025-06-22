# =============================================================================
# TileMiner — Sequence Tile Plotter for NCBI Nucleotide Data
# -----------------------------------------------------------------------------
# Author: Jason Fajardo
#
# Description:
#     TileMiner is a customizable R tool for fetching and visualizing nucleotide 
#     sequences from NCBI as checker-style tile plots. This helps reveal 
#     layout-dependent structure and pattern emergence across biological sequences.
#
# Features:
#     - Supports single accessions and Entrez batch queries
#     - Checker-style visual tile plots with user-defined column widths
#     - Sequence slicing support for accessions
#     - User-defined color schemes for tile rendering
#     - Metadata extraction (e.g., gene name, organism, length)
#     - Optional FASTA export and tile pagination
#     - Automatic plot export for all widths and pages
#     - Summary reports and per-sequence metadata files (.meta)
#
# Requirements:
#     R packages: rentrez, Biostrings, ggplot2, reshape2
#
# License:
#     Code: MIT License
#     TileSuit Visual Language: Creative Commons CC BY 4.0
#
# -----------------------------------------------------------------------------
# === Usage Overview ===
#
# tile_miner("NM_001101.5")  
# tile_miner("actin[gene] AND Homo sapiens[Organism]")  
#
# Output:
#   - PNG plots for each sequence (one per width, per page)
#   - Optional: FASTA file per sequence
#   - Per-sequence metadata (.meta)
#   - Summary files: summary.txt and summary.csv
#
# -----------------------------------------------------------------------------
# === tile_miner() — Universal Wrapper Function ===
#
# Automatically detects whether input is:
#   - A direct accession (e.g., "NM_001101.5") → calls tile_fetch_one()
#   - An Entrez query (e.g., "actin[Gene] AND Homo sapiens[Organism]") → calls tile_fetch_batch()
#
# Parameters:
#   query           — Accession or NCBI Entrez query string
#   widths          — Vector of tile column widths (e.g., c(40, 60))
#   bases_per_page  — Number of bases per tile plot page (default: 5000)
#   min_length      — Filter: minimum sequence length
#   max_length      — Filter: maximum sequence length
#   save_fasta      — Export sequences as FASTA (TRUE/FALSE)
#   verbose         — Print progress output
#   sort            — Sort method for Entrez query results (e.g., "relevance", "date")
#   max_results     — Maximum number of sequences to retrieve (batch mode)
#   batch_size      — Number of sequences processed per batch
#   start, end      — (Optional) Start and end coordinates for slicing (accession mode only)
#   color_scheme    — (Optional) Named vector or function for overriding default color mapping
#
# Notes:
#   - Coordinates are 1-based and inclusive (like Biostrings).
#   - Layout (`widths`) greatly affects visual structure—resonant patterns can emerge.
#   - Custom color schemes allow domain-specific encodings (e.g., GC skew, entropy)
#
# -----------------------------------------------------------------------------
# === Example Calls ===
#
# Visualize full sequence:
#   tile_miner("NM_001101.5", widths = c(60))
#
# Slice a genome region:
#   tile_miner("NC_000913.3", start = 363231, end = 366305, widths = 60)
#
# Batch query with filters:
#   tile_miner("actin[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
#              max_results = 5)
#
# -----------------------------------------------------------------------------
# === Extended Test Cases ===
#
# 1. Multiple widths, large sequence:
#    tile_miner("NC_045512.2", widths = c(40, 60, 80), bases_per_page = 2000)
#
# 2. FASTA disabled, length filtered:
#    tile_miner("NM_001101.5", widths = 50, save_fasta = FALSE, 
#               min_length = 1000, max_length = 3000)
#
# 3. Batch query, small batch size:
#    tile_miner("insulin[Title] AND biomol_mrna[PROP]", max_results = 3,
#               batch_size = 1, widths = c(40, 80))
#
# 4. Filter by publication date:
#    tile_miner("p53[Gene] AND Homo sapiens[Organism] AND srcdb_refseq[PROP] AND 2010/01/01:2024/12/31[PDAT]",
#               widths = c(60), max_results = 5, min_length = 1000)
#
# 5. Graceful failure (overly strict):
#    tile_miner("unknown_gene_name[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
#               max_results = 5)
#
# 6. Bacterial genomic gene:
#    tile_miner("rpoB[Gene] AND Escherichia coli[Organism] AND biomol_genomic[PROP]",
#               widths = 70, max_results = 3, min_length = 5000)
#
# 7. RNA virus genome:
#    tile_miner("NC_001802.1", widths = c(40, 60, 100), bases_per_page = 1500)
#
# -----------------------------------------------------------------------------
# === NCBI Filter Cheatsheet (for queries) ===
#
#   [Organism]         => Homo sapiens[Organism]
#   [Gene]             => actin[Gene], p53[Gene]
#   [Title]            => Less strict alternative to [Gene]
#   [Sequence Length]  => 100:2000[Sequence Length]
#   [PDAT]             => Publication date, e.g., 2020/01/01:2024/12/31[PDAT]
#
#   [PROP] Filters:
#     biomol_mrna[PROP]     → mRNA
#     biomol_genomic[PROP]  → Genomic DNA
#     biomol_rrna[PROP]     → rRNA
#     srcdb_refseq[PROP]    → RefSeq source
#
# Example Query Combos:
#   - "actin[Title] AND Homo sapiens[Organism] AND biomol_mrna[PROP]"
#   - "p53[Gene] AND Mus musculus[Organism] AND srcdb_refseq[PROP]"
#   - "insulin[Title] AND 1000:3000[Sequence Length]"
#
# Tip: Start broad, then tighten filters. Over-constraining will often yield zero results.
# =============================================================================

#
#
# Required packages
library(rentrez)
library(Biostrings)
library(ggplot2)
library(reshape2)

# =====================================================================================================
# === BLOCK: Helper Functions ===

# These internal functions support the tile plotting pipeline.

# You can customize and tweak behavior by adjusting them as noted below.

# Checks if input is a valid NCBI-style accession number (e.g., NM_001101.5)
is_accession <- function(x) {
  grepl("^[A-Z]{1,2}_[0-9]+(\\.[0-9]+)?$", x) && !grepl(" ", x)
}

# Removes all non-nucleotide characters and capitalizes valid bases
# You can tweak the character filter if using sequences with ambiguity codes
clean_seq <- function(seq) {
  toupper(gsub("[^ATGCUN]", "", seq))
}

# Writes the sequence to a FASTA file using accession as the filename
# The header contains accession + gene + full title
save_fasta_file <- function(accession, gene, full_name, seq, outdir) {
  label <- accession
  header <- paste0(">", accession, if (!is.na(gene) && nzchar(gene)) paste0(" ", gene) else "", " ", full_name)
  fasta_text <- paste0(header, "\n", seq)
  writeLines(fasta_text, file.path(outdir, paste0(accession, ".fasta")))
}

# Converts a sequence character vector into a tile-ready data frame
# Customize row/column offset behavior here if changing visual layout
# Ensure 'width' matches plot layout column count
to_matrix_df <- function(chars, width) {
  len <- length(chars)
  height <- ceiling(len / width)
  chars <- c(chars, rep(NA, height * width - len))
  mat <- matrix(chars, nrow = height, byrow = TRUE)
  df <- melt(mat, varnames = c("row", "col"), value.name = "char")
  df$row_offset <- max(df$row) - df$row + 1
  df$col_offset <- df$col
  df <- df[!is.na(df$char), ]
  attr(df, "tile_matrix") <- mat
  return(df)
}

# =====================================================================================================

# === render_tile_plot(): Renders a single tile plot using ggplot2. ===

# You can adjust:
# - Title format: lines 2–3 below
# - Color scheme: inside scale_fill_manual()
# - Tile size/layout: coord_fixed() and theme()
# - Output size: ggsave(... width = , height = )

# colour schemes can be preset by creating a vector of colours under color_scheme
# i.e.,
#color_scheme <- c(A = "green", T = "red", C = "blue", G = "pink", N = "orange")

# and passing the function as an argument for tile_miner() or any of the other functions.

# === BLOCK: render_tile_plot() ===
render_tile_plot <- function(df, width, page, accession, outdir, seq_len, total_pages, full_name = NULL, color_scheme = NULL) {
  if (is.null(color_scheme)) {
    color_scheme <- c(A = "black", T = "gray30", C = "gray30", G = "black", N = "orange")
  }
  title <- if (!is.null(full_name)) paste0(accession, "\n", full_name) else accession
  subtitle <- paste("Width:", width, "| Page:", page, "of", total_pages)
  caption <- paste0("Length: ", seq_len, "bp\nColumns: ", width, "\nTotal Pages: ", total_pages)
  g <- ggplot(df, aes(x = col_offset, y = row_offset)) +
    geom_tile(aes(fill = char), colour = "white", linewidth = 0) +
    scale_fill_manual(values = color_scheme, na.value = "white") +
    theme_void() +
    coord_fixed() +
    ggtitle(title, subtitle = subtitle) +
    labs(caption = caption) +
    theme(plot.caption = element_text(hjust = 0, size = 8),
          plot.subtitle = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold"))
  
  # Save PNG
  png_filename <- paste0(gsub("[[:punct:]]", "_", accession), "_w", width, "_p", page, ".png")
  ggsave(file.path(outdir, png_filename), g, width = 10, height = 10)
  
  # Save metadata as .meta
  meta_filename <- paste0(gsub("[[:punct:]]", "_", accession), "_w", width, "_p", page, ".meta")
  
  meta <- list(
    accession = accession,
    width = width,
    page = page,
    sequence_length = seq_len,
    layout = list(
      rows = max(df$row_offset),
      cols = max(df$col_offset)
    ),
    matrix = attr(df, "tile_matrix")
  )
  saveRDS(meta, file = file.path(outdir, meta_filename))
}



# === BLOCK: validate_meta_files() ===
# Utility function to test whether all .rds metadata files in a folder can be successfully read
validate_meta_files <- function(folder) {
  rds_files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("No .rds files found in: ", folder)
    return(data.frame(
      file = character(0),
      readable = logical(0),
      error = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  results <- data.frame(
    file = basename(rds_files),
    readable = logical(length(rds_files)),
    error = character(length(rds_files)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(rds_files)) {
    tryCatch({
      readRDS(rds_files[i])
      results$readable[i] <- TRUE
    }, error = function(e) {
      results$readable[i] <- FALSE
      results$error[i] <- conditionMessage(e)
    })
  }
  
  return(results)
}


# =====================================================================================================

# === process_accession(): Core function for fetching, parsing, and plotting one accession ===

# Handles all internal logic: metadata extraction, sequence cleaning, pagination, tile rendering.
# You can customize:
# - How sequences are sliced per page (bases_per_page)
# - Visual layout options (widths)
# - Plot file naming/formatting

process_accession <- function(id, outdir, widths, bases_per_page, min_length, max_length, save_fasta, verbose, color_scheme = NULL) {
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
  if (save_fasta) save_fasta_file(accession, gene, full_name, seq, outdir)
  chars <- strsplit(seq, "")[[1]]
  for (w in widths) {
    total_len <- length(chars)
    
    # Adjust bases_per_page to ensure full rows per page by flooring the division
    # This guarantees consistent tile height and avoids visual truncation
    
    rows_per_page <- floor(bases_per_page / w)
    bases_per_page_adj <- rows_per_page * w
    pages <- ceiling(total_len / bases_per_page_adj)
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
      start <- ((p - 1) * bases_per_page_adj) + 1
      end <- min(start + bases_per_page_adj - 1, total_len)
      df <- to_matrix_df(chars[start:end], w)
      render_tile_plot(df, w, p, accession, outdir, total_len, pages, full_name = header_metadata, color_scheme = color_scheme)
    }
  }
  return(metadata_summary)
}

# =====================================================================================================

# === BLOCK: tile_fetch_one(): Process a single accession ===

# Entry point for direct plotting of one NCBI accession.
# You can tweak:
# - widths: Number of columns in tile plots (affects layout/resolution)
# - bases_per_page: Controls how much sequence per page (tune for pattern visibility)
# - min_length/max_length: Filter sequences by size before plotting
# - save_fasta: Whether to store FASTA files
# - verbose: Print progress messages during run

tile_fetch_one <- function(accession,
                           widths = c(60),
                           bases_per_page = 5000,
                           min_length = 200,
                           max_length = NULL,
                           save_fasta = TRUE,
                           verbose = TRUE,
                           start = NULL,
                           end = NULL,
                           color_scheme = NULL) {
  if (verbose) message("Detected accession input — skipping search.")
  outdir <- paste0(gsub("[[:punct:]]", "_", accession), "_tileminer")
  dir.create(outdir, showWarnings = FALSE)
  
  # Fetch full sequence
  fasta <- entrez_fetch(db = "nuccore", id = accession, rettype = "fasta", retmode = "text")
  lines <- strsplit(fasta, "\\n")[[1]]
  header <- gsub("^>", "", lines[1])
  seq_full <- clean_seq(paste(lines[-1], collapse = ""))
  
  # If slicing coordinates provided, bypass process_accession and handle inline
  if (!is.null(start) && !is.null(end)) {
    seq_sliced <- substr(seq_full, start, end)
    chars <- strsplit(seq_sliced, "")[[1]]
    sliced_accession <- paste0(accession, "_", start, "-", end)
    full_name <- paste(header)
    metadata_summary <- list()
    csv_rows <- list()
    for (w in widths) {
      total_len <- length(chars)
      rows_per_page <- floor(bases_per_page / w)  # Floors rows for consistent full tile rows
      bases_per_page_adj <- rows_per_page * w     # Adjusts actual bases per page
      pages <- ceiling(total_len / bases_per_page_adj)
      entry <- list(
        accession = sliced_accession,
        name = full_name,
        length = total_len,
        width = w,
        pages = pages
      )
      metadata_summary[[length(metadata_summary) + 1]] <- entry
      for (p in seq_len(pages)) {
        if (verbose) message("Plotting ", sliced_accession, " | Width ", w, " | Page ", p)
        start_idx <- ((p - 1) * bases_per_page_adj) + 1
        end_idx <- min(start_idx + bases_per_page_adj - 1, total_len)
        df <- to_matrix_df(chars[start_idx:end_idx], w)
        render_tile_plot(df, w, p, sliced_accession, outdir, total_len, pages, full_name = full_name, color_scheme = color_scheme)
      }
    }
  } else {
    metadata_summary <- process_accession(accession, outdir, widths, bases_per_page, min_length, max_length, save_fasta, verbose, color_scheme = color_scheme)
  }
  
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

# =====================================================================================================

# === BLOCK: tile_fetch_batch(): Perform search and batch fetch ===

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

# === tile_fetch_batch(): Perform search and batch fetch ===

# Entry point for Entrez search + multi-sequence plotting.
# You can customize:
# - sort: Order of returned search hits (e.g., by length, date)
# - max_results: Max number of hits to retrieve
# - batch_size: Controls how many are processed at a time (avoid timeouts)
# - widths, bases_per_page, min/max length, etc. same as tile_fetch_one

tile_fetch_batch <- function(query,
                             sort = "relevance",
                             widths = c(60),
                             max_results = 10,
                             batch_size = 5,
                             bases_per_page = 5000,
                             min_length = 200,
                             max_length = NULL,
                             save_fasta = TRUE,
                             verbose = TRUE,
                             color_scheme = NULL) {
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
        result <- process_accession(id, outdir, widths, bases_per_page, min_length, max_length, save_fasta, verbose, color_scheme = color_scheme)
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


# === BLOCK: Wrapper Function ===
# - Supports all parameters from tile_fetch_one() and tile_fetch_batch(), including:
#     widths, max_results, batch_size, bases_per_page, min_length, max_length,
#     save_fasta, verbose, and sort (controls NCBI search order).
#     Valid sort values: "relevance" (default), "date", "submit_date", "sequence_length", "accn", "source"
# === tile_miner(): Universal wrapper function ===
# Automatically detects if query is a single accession or a search term.
# For accessions, calls tile_fetch_one(). For queries, calls tile_fetch_batch().
# All parameters (widths, batch_size, etc.) are passed through.
tile_miner <- function(query,
                       widths = c(60),
                       bases_per_page = 5000,
                       min_length = 200,
                       max_length = NULL,
                       save_fasta = TRUE,
                       verbose = TRUE,
                       sort = "relevance",
                       max_results = 10,
                       batch_size = 5,
                       start = NULL,
                       end = NULL,
                       color_scheme = NULL) {
  if (is_accession(query)) {
    tile_fetch_one(query,
                   widths = widths,
                   bases_per_page = bases_per_page,
                   min_length = min_length,
                   max_length = max_length,
                   save_fasta = save_fasta,
                   verbose = verbose,
                   start = start,
                   end = end,
                   color_scheme = color_scheme)
  } else {
    tile_fetch_batch(query,
                     sort = sort,
                     widths = widths,
                     max_results = max_results,
                     batch_size = batch_size,
                     bases_per_page = bases_per_page,
                     min_length = min_length,
                     max_length = max_length,
                     save_fasta = save_fasta,
                     verbose = verbose,
                     color_scheme = color_scheme)
  }
}

# =====================================================================================================

# === Extended Test Examples ===

# These test a wide range of parameters and scenarios:

# 1. Single accession, large sequence, multiple widths: (VERY LARGE)
# tile_miner("NC_045512.2", widths = c(40, 60, 80), bases_per_page = 2000)

# 2. Single accession, with FASTA disabled and custom length limits:
# tile_miner("NM_001101.5", widths = 50, save_fasta = TRUE, min_length = 1000, max_length = 3000)

# 3. Batch query, max 3 results, small batch size, two widths:
# tile_miner("insulin[Title] AND biomol_mrna[PROP]", max_results = 3, batch_size = 1, widths = c(40, 80))

# 4. Longer query with organism, source, and publication date filter:
#tile_miner("p53[Gene] AND Homo sapiens[Organism] AND srcdb_refseq[PROP] AND 2010/01/01:2024/12/31[PDAT]",
  # widths = c(60), max_results = 5, min_length = 1000)

# 5. Intentionally too-strict query (tests graceful failure):
# tile_miner("unknown_gene_name[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
  # max_results = 5)

# 6. Bacterial gene, medium size, with genomic source filter:
# tile_miner("rpoB[Gene] AND Escherichia coli[Organism] AND biomol_genomic[PROP]",
  # widths = 70, max_results = 50, min_length = 5000)

# 7. RNA viral sequence, small genome:
# tile_miner("NC_001802.1")

# 8. Mitochondrial sequence.
# tile_fetch_batch("mitochondrion[Filter] AND Homo sapiens[Organism]")

# 9. Multiple Sequences One Call, Three Plots.

tile_miner(
  query = "NM_001101.5",
  widths = c(60, 80, 90, 100),
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


# 11. Slice a specific region from a long genome:
tile_miner("NC_000913.3", start = 363231, end = 366305, widths = 60, bases_per_page = 3000)


# === Examples to Try Colour Scheme ===

# 1. tile_miner() — top-level wrapper
tile_miner("NM_001101.5", widths = 50, color_scheme = c(A = "red", T = "blue", G = "green", C = "yellow"))

# 2. tile_fetch_one() — direct accession with slicing
# Use case: Slice small region and apply high-contrast scheme

tile_fetch_one(
  "NM_007393.5",
  start = 100,
  end = 400,
  widths = 30,
  color_scheme = c(A = "#ffaaaa", T = "#aaaaff", G = "#aaffaa", C = "#ffd700")
)

# 3. tile_fetch_batch() — search with multiple small accessions
# Use case: Batch with muted grayscale palette for entropy analysis

tile_fetch_batch(
  query = "insulin[Title] AND Homo sapiens[Organism] AND biomol_mrna[PROP]",
  widths = 50,
  max_results = 3,
  batch_size = 1,
  color_scheme = c(A = "gray20", T = "gray40", G = "gray60", C = "gray80")
)



# =======================================================================================================
# Meta file validation
validate_meta_files("actin_beta_Title__AND__Homo_sapiens_Orga_tileminer")

readRDS("NM_001101_5_tileminer/NM_001101_5_w50_p1.meta")
readRDS("mitochondrion_Filter__AND_Homo_sapiens_O_tileminer/PV695974_1_w60_p2.meta")

