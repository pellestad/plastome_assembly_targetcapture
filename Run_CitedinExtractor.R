setwd("~/path/to/your/project")

# *** MUST be the very first thing that runs, in a completely fresh R session ***

Sys.setenv(CURL_SSL_BACKEND = "openssl")

options(warn = 1)

# Load packages
library(rentrez)
library(easyPubMed)
library(httr)

# Load the CitedInExtractoR() function
# (make sure CitedInExtractoR.R is in your working directory, or give its full path)
source("CitedInExtractoR.R")

# reduces how often a connection gets reused after the server half-closed it.
httr::set_config(httr::config(forbid_reuse = 1L, http_version = 1.1))

# Optional but recommended: register an NCBI API key to raise the rate limit
# from 3 to 10 requests/second (get one at https://www.ncbi.nlm.nih.gov/account/settings/).
ncbi_key <- Sys.getenv("NCBI_API_KEY")
if (nzchar(ncbi_key)) {
  rentrez::set_entrez_key(ncbi_key)
} else {
  message("No NCBI_API_KEY environment variable set -- continuing at the ",
          "unauthenticated 3 requests/sec rate limit.")
}

# List PubMed IDs for the target publications you want citation data for
PubMedID <- c("32912315", "28204566", "23661685")

# Apply function across PubMedID vector
DatPubMed <- NULL
for(i in seq_along(PubMedID)){
  tmp <- tryCatch(
    CitedInExtractoR(PubMedID[i]),
    error = function(e){
      warning(paste("CitedInExtractoR failed for", PubMedID[i], ":", conditionMessage(e)))
      NULL
    }
  )
  if(!is.null(tmp)) DatPubMed <- rbind(DatPubMed, tmp)
}

# Sanity check before assuming column names: table_articles_byAuth()'s exact
# output columns have drifted across easyPubMed versions/docs, so confirm
# what's actually here before relying on a "doi" column.
if(!is.null(DatPubMed)){
  print(names(DatPubMed))
  print(head(DatPubMed))
}

# Add a URL to directly access each DOI (added once, after the full table is built)
# - only if a doi-like column actually exists under one of its known names.
if(!is.null(DatPubMed)){
  doi_col <- intersect(c("doi", "article.doi", "pmid.doi"), names(DatPubMed))[1]
  if(!is.na(doi_col)){
    DatPubMed$URLDOI <- paste0("http://doi.org/", DatPubMed[[doi_col]])
  } else {
    warning("No doi-like column found - check the printed column names above and adjust doi_col.")
  }
}

# Write data out (as csv)
write.csv(
  DatPubMed,
  file = "PubCitedIn_publications.csv",
  row.names = FALSE,
  quote = TRUE
)
