#######################################
#    Plastome Seed Finder Pipeline    #
#######################################
#   By John M. A. Wojahn, PhD FLS     #
#######################################
# Released under the CC BY-NC-BY 4.0  #
#######################################
#         April 24th, 2026            #
#######################################
## Reference plastome selection
#
#Picks a "seed" reference plastome for each sample in a sample sheet, for use
#as the reference genome in reference-guided steps
#
#For each sample it tries, in order:
#  
#1. An exact species-level search against NCBI `nuccore` for a complete, annotated plastid/chloroplast genome.
#2. If nothing at the species level, a genus-level search.
#3. If nothing at the genus level either, a phylogenetic fallback: starting from the sample's genus tip on a PAFTOL reference tree, 
#it walks outward through successive sister clades and tries each sister taxon name against NCBI in turn, stopping at the first one 
#with a valid hit (or after visiting 200 tips, or if the genus isn't in the tree at all).
#
# ---------------------------------------------------------------------------
# Required inputs
# ---------------------------------------------------------------------------
#
# - A PAFTOL reference tree in Newick format (--tree; default
#   `treeoflife.current.tree`), whose tip labels are underscore-delimited
#   with the genus name as the third field (e.g.
#   `something_something_Genus_something`) -- the script strips them down
#   to just that genus field on load, same as the original.
# - Your sample sheet (--sampling; see format above).
#
# Requirements: R packages `rentrez`, `ape`, `xml2`.
# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
#
#   Rscript find_reference_plastomes.R \
#       --tree treeoflife.current.tree \
#       --sampling sampling.csv \
#       --genus-col genus \
#       --species-col species \
#       --min-length 100000 \
#       --out samplez_done_new.csv
#
# All flags are optional; the values above are also the defaults, so for a
# sampling.csv that already has "genus" and "species" columns, just:
#
#   Rscript find_reference_plastomes.R
#
# ---------------------------------------------------------------------------
# Expected sampling.csv format
# ---------------------------------------------------------------------------
#
# A CSV with (at minimum) one column giving each sample's genus and one
# giving a species-level search term, e.g.:
#
#   sample_id,genus,species,notes
#   S001,Bidens,Bidens alba,
#   S002,Calendula,Calendula arvensis,
#
# - "genus" (or whatever column you point --genus-col at) must contain just
#   the genus name, e.g. "Bidens" -- it's matched directly against the
#   (genus-only) tip labels of the PAFTOL tree for the phylogenetic
#   fallback, so it needs to use the same genus spelling/synonymy as the
#   tree.
# - "species" (or --species-col) is used as-is as an NCBI nuccore
#   [Organism] search term for the first (species-level) search attempt --
#   typically a full binomial like "Bidens alba", but anything that's a
#   valid NCBI organism-field search term works.
# - Any other columns (sample IDs, notes, etc.) are carried through
#   untouched into the output CSV.
# - Row order and row count are not assumed -- every row in the file is
#   processed, and rows are no longer silently deduplicated or truncated
#   (the original script's `samplez[1:132, 1:4]` slice and blanket
#   `unique()` call, both specific to the original 132-row sheet, have been
#   removed -- dedupe your own input first if you want that).
#
# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
#Results (including `MANUAL_NOT_IN_TREE` / `MANUAL_NOT_FOUND` / `ERROR` placeholders for samples nothing was found for) are written back out
#alongside the original sample sheet columns, prefixed `ncbi_` (see "Output" below).


suppressMessages(library(rentrez))
suppressMessages(library(ape))
suppressMessages(library(xml2))

# ---------------------------------------------------------------------------
# CLI ARGUMENT PARSING
# ---------------------------------------------------------------------------
# Simple `--flag value` parser -- no extra package dependency (e.g.
# optparse) beyond what this script already required.
parse_cli_args <- function(args, defaults) {
  cfg <- defaults
  flag_to_key <- c(
    "--tree"        = "tree",
    "--sampling"    = "sampling",
    "--genus-col"   = "genus_col",
    "--species-col" = "species_col",
    "--min-length"  = "min_length",
    "--out"         = "out"
  )
  i <- 1
  while (i <= length(args)) {
    flag <- args[i]
    if (!(flag %in% names(flag_to_key))) {
      stop("Unrecognized argument: ", flag,
           "\nRecognized flags: ", paste(names(flag_to_key), collapse = ", "))
    }
    key <- flag_to_key[[flag]]
    if (i == length(args)) {
      stop("Flag ", flag, " is missing its value.")
    }
    value <- args[i + 1]
    if (key == "min_length") value <- suppressWarnings(as.numeric(value))
    cfg[[key]] <- value
    i <- i + 2
  }
  cfg
}

defaults <- list(
  tree        = "treeoflife.current.tree",
  sampling    = "sampling.csv",
  genus_col   = "genus",
  species_col = "species",
  min_length  = 100000, # Min threshold; 100kb-80kb is a reasonable range for a plastome
  out         = "samplez_done_new.csv"
)
CFG <- parse_cli_args(commandArgs(trailingOnly = TRUE), defaults)

if (is.na(CFG$min_length)) {
  stop("--min-length must be a number (got: ", defaults$min_length, ")")
}
if (!file.exists(CFG$tree)) {
  stop("Tree file not found: ", CFG$tree,
       "\n(point --tree at your PAFTOL Newick tree)")
}
if (!file.exists(CFG$sampling)) {
  stop("Sampling sheet not found: ", CFG$sampling,
       "\n(point --sampling at your sample CSV)")
}

# ---------------------------------------------------------------------------
# READ INPUTS
# ---------------------------------------------------------------------------
paftol_tree <- ape::read.tree(CFG$tree)
newlab <- rep(NA, length(paftol_tree$tip.label))
for (i in seq_along(paftol_tree$tip.label)) {
  newlab[i] <- unlist(strsplit(paftol_tree$tip.label[i], split = "_"))[3]
}
paftol_tree$tip.label <- newlab

samplez <- read.csv(CFG$sampling, stringsAsFactors = FALSE)

missing_cols <- setdiff(c(CFG$genus_col, CFG$species_col), names(samplez))
if (length(missing_cols) > 0) {
  stop("Column(s) not found in ", CFG$sampling, ": ",
       paste(missing_cols, collapse = ", "),
       "\nColumns present: ", paste(names(samplez), collapse = ", "),
       "\n(use --genus-col / --species-col to point at the right columns)")
}

# NCBI API key: an earlier version of this script had a real key hardcoded
# directly in the source (`rentrez::set_entrez_key("...")`). That key has
# been removed before adding this script to the public repo -- if you're
# the owner of that key, consider it exposed and generate a fresh one from
# your NCBI account settings (https://www.ncbi.nlm.nih.gov/account/settings/).
# Registering a free key raises the E-utilities rate limit from 3 to 10
# requests/second. Don't hardcode your own key here either -- set it via the
# NCBI_API_KEY environment variable instead (e.g. in your ~/.Renviron).
ncbi_key <- Sys.getenv("NCBI_API_KEY")
if (nzchar(ncbi_key)) {
  rentrez::set_entrez_key(ncbi_key)
} else {
  message("No NCBI_API_KEY environment variable set -- continuing at the ",
          "unauthenticated 3 requests/sec rate limit.")
}

uwu_factory <- function(paftol_tree, min_length = CFG$min_length) {

  # --- INTERNAL UTILS ---
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # --- CACHES ---
  seq_cache <- new.env(parent = emptyenv())
  tax_cache <- new.env(parent = emptyenv())

  # --- HELPER 1: CHUNKED NCBI FETCH (Fixes Error 414) ---
  get_seq_summaries <- function(ids) {
    ids <- unique(ids)
    cached <- ids[ids %in% ls(seq_cache)]
    new_ids <- setdiff(ids, cached)
    out <- list()
    if (length(cached) > 0) out[cached] <- mget(cached, envir = seq_cache)

    if (length(new_ids) > 0) {
      chunk_size <- 40 # Conservative chunk size for long URLs
      id_chunks <- split(new_ids, ceiling(seq_along(new_ids) / chunk_size))
      for (chunk in id_chunks) {
        sums <- entrez_summary(db = "nuccore", id = chunk)
        if (length(chunk) == 1) {
            sums <- setNames(list(sums), chunk)
        }
        for (id in names(sums)) {
          seq_cache[[id]] <- sums[[id]]
          out[[id]] <- sums[[id]]
        }
        Sys.sleep(0.3)
      }
    }
    return(out)
  }

  # --- HELPER 2: SEARCHER (The missing function) ---
  basicsearcheR <- function(term) {
    term <- trimws(term)
    if (is.na(term) || term == "") return(list(count = 0, ids = character(0)))
    res <- tryCatch({
      rentrez::entrez_search(
        db = "nuccore",
        term = paste0(term, "[Organism] AND (biomol_genomic[PROP] AND complete[TITLE])"),
        use_history = TRUE,
        retmax = 500
      )
    }, error = function(e) {
      message("Search failed for ", term, ": ", e$message)
      return(list(count = 0, ids = character(0)))
    })
    Sys.sleep(0.3)
    list(count = res$count %||% 0, ids = as.character(res$ids %||% character(0)))
  }

  # --- HELPER 3: TAXONOMY ---
  get_taxonomy <- function(taxid) {
    if (is.null(taxid) || is.na(taxid)) return(NULL)
    taxid <- as.character(taxid)
    if (exists(taxid, envir = tax_cache)) return(tax_cache[[taxid]])
    xml <- tryCatch({
      entrez_fetch(db = "taxonomy", id = taxid, rettype = "xml", parsed = FALSE)
    }, error = function(e) NULL)
    if (is.null(xml)) return(NULL)
    xml_parsed <- xml2::read_xml(xml)
    tax_cache[[taxid]] <- xml_parsed
    return(xml_parsed)
  }

  # --- HELPER 4: PARSE TAX ---
  parse_tax <- function(xml) {
    ranks <- c("kingdom","phylum","class","order","family","genus","species")
    out <- setNames(as.list(rep(NA, length(ranks))), ranks)
    if (is.null(xml)) return(out)
    nodes <- xml2::xml_find_all(xml, ".//Taxon")
    for (n in nodes) {
      rank <- xml2::xml_text(xml2::xml_find_first(n, "./Rank"))
      name <- xml2::xml_text(xml2::xml_find_first(n, "./ScientificName"))
      if (rank %in% ranks) out[[rank]] <- name
    }
    sp <- xml2::xml_text(xml2::xml_find_first(xml, ".//ScientificName"))
    if (!is.na(sp) && sp != "") out$species <- sp
    return(out)
  }

  # --- HELPER 5: VALIDATE & ENRICH ---
  get_valid_records <- function(ids) {
    if (length(ids) == 0) return(NULL)
    seqs <- get_seq_summaries(ids)
    res_list <- lapply(names(seqs), function(id) {
      s <- seqs[[id]]
      if (is.null(s)) return(NULL)
      slen <- s$slen %||% s$length %||% NA
      genome_type <- s$genome %||% NA
      comp <- s$completeness %||% NA
      if (!is.na(genome_type) &&
          genome_type %in% c("chloroplast","plastome","plastid") &&
          !is.na(comp) && comp == "complete" &&
          !is.na(slen) && slen >= min_length) {
        tax_data <- parse_tax(get_taxonomy(as.character(s$taxid)))
        return(cbind(data.frame(id=id, genome=genome_type, length=slen, taxid=as.character(s$taxid)), as.data.frame(tax_data)))
      }
      return(NULL)
    })
    # NOTE: the closing `})` above (for the lapply() call started a few
    # lines up) was missing in the original script -- as written, the
    # anonymous function and the lapply() call were never closed, which is
    # a guaranteed R parse error before this function could ever run. Added
    # to fix that; nothing about the matching logic itself was changed.
    df <- do.call(rbind, res_list)
    return(df)
  }

  # --- HELPER 6: PHYLO SISTER LOCATOR ---
  get_sister_tips <- function(tree, node) {
    edge_row <- which(tree$edge[, 2] == node)
    if (length(edge_row) == 0) return(character(0))
    parent <- tree$edge[edge_row, 1]
    sisters <- tree$edge[tree$edge[, 1] == parent, 2]
    sisters <- setdiff(sisters, node)
    tips <- unlist(lapply(sisters, function(x) {
      if (x <= length(tree$tip.label)) tree$tip.label[x]
      else ape::extract.clade(tree, x)$tip.label
    }))
    return(unique(tips))
  }

  # ==========================================
  # THE RETURNED ANALYZER FUNCTION
  # ==========================================
  # Takes the genus and species search term for one sample directly (rather
  # than a whole sample-sheet + row index and a hardcoded column position),
  # so it works regardless of what your sampling.csv's other columns look
  # like.
  function(genus, species_term) {
    message(sprintf("Processing: %s (genus: %s)", species_term, genus))
    tryCatch({
      # 1. Species Match
      res <- basicsearcheR(species_term)
      if (res$count > 0) {
        df <- get_valid_records(res$ids)
        if (!is.null(df)) return(as.list(df[sample.int(nrow(df), 1), ]))
      }
      # 2. Genus Match
      res <- basicsearcheR(genus)
      if (res$count > 0) {
        df <- get_valid_records(res$ids)
        if (!is.null(df)) return(as.list(df[sample.int(nrow(df), 1), ]))
      }
      # 3. Phylo Fallback
      visited_tips <- genus
      current_node <- which(paftol_tree$tip.label == genus)
      if (length(current_node) == 0) return(list(id="MANUAL_NOT_IN_TREE"))

      repeat {
        sisters <- get_sister_tips(paftol_tree, current_node)
        sisters <- setdiff(sisters, visited_tips)
        if (length(sisters) > 0) {
          for (cand in sample(sisters)) {
            res <- basicsearcheR(cand)
            if (res$count == 0) next
            df <- get_valid_records(res$ids)
            if (!is.null(df)) return(as.list(df[sample.int(nrow(df), 1), ]))
            visited_tips <- c(visited_tips, cand)
          }
        }
        parent_row <- which(paftol_tree$edge[, 2] == current_node)
        if (length(parent_row) == 0) break
        current_node <- paftol_tree$edge[parent_row, 1]
        if (length(visited_tips) > 200) break
      }
      return(list(id="MANUAL_NOT_FOUND"))
    }, error = function(e) {
      message("Error processing ", species_term, ": ", e$message)
      return(list(id="ERROR"))
    })
  }
}

# ---------------------------------------------------------------------------
# EXECUTION
# ---------------------------------------------------------------------------
AnalylizeR <- uwu_factory(paftol_tree, CFG$min_length)
results <- vector("list", nrow(samplez))
for (i in seq_len(nrow(samplez))) {
  results[[i]] <- AnalylizeR(samplez[[CFG$genus_col]][i], samplez[[CFG$species_col]][i])
}

# Standardize column names across the list before rbinding
# (Some might have 'id' only if they failed)
all_cols <- c("id", "genome", "length", "taxid", "kingdom", "phylum", "class", "order", "family", "genus", "species")
results_cleaned <- lapply(results, function(x) {
  for (col in all_cols) if (!col %in% names(x)) x[[col]] <- NA
  return(as.data.frame(x))
})

results_df <- do.call(rbind, results_cleaned)
# Prefix the matched-record columns before merging them onto your sample
# sheet. Without this, a sampling.csv that (very plausibly, for a taxonomic
# sample sheet) already has its own "genus", "species", or "family" columns
# would end up with duplicate column names after cbind() -- confusing, and
# liable to silently grab the wrong one if you later do samplez_done$genus.
names(results_df) <- paste0("ncbi_", names(results_df))
samplez_done <- cbind(samplez, results_df)
write.csv(samplez_done, CFG$out, row.names = FALSE)
message("Wrote ", nrow(samplez_done), " row(s) to ", CFG$out)
