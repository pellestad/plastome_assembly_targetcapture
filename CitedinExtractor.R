# Generic retry wrapper. NCBI's E-utilities occasionally drop the TLS
# connection mid-response - on Windows this surfaces via schannel as
# "server closed abruptly (missing close_notify)". It's a transient/benign
# hiccup (curl's strict TLS shutdown check, not lost data), so we just retry
# with increasing backoff instead of letting it kill the whole run.
.retry <- function(fn, max_retries = 4, delay = 1, label = "request"){
  attempt <- 1
  repeat{
    result <- tryCatch(list(ok = TRUE, value = fn()), error = function(e){
      list(ok = FALSE, err = conditionMessage(e))
    })
    if(result$ok) return(result$value)
    if(attempt >= max_retries){
      warning(paste(label, "failed after", attempt, "attempts:", result$err))
      return(NULL)
    }
    warning(paste(label, "attempt", attempt, "failed (", result$err, "), retrying..."))
    Sys.sleep(delay * attempt)
    attempt <- attempt + 1
  }
}

#A function to list publications citing a target article and
# downloading data on these latter publications
CitedInExtractoR <- function(PubMedID, batch_size = 50, delay = 0.4, max_retries = 4){
  print(paste("Fetch data associated to", PubMedID, sep = " "))
  #Search data associated to PubMedID in pubmed (retried - this call can drop mid-request)
  # linkname restricts elink to computing ONLY the citedin links instead of every
  # neighbor category (related articles, PMC refs, etc.) - a much smaller/faster
  # response, which is less likely to get cut off mid-transfer.
  src <- .retry(
    function() rentrez::entrez_link(dbfrom = "pubmed", id = PubMedID, db = "pubmed",
                                    linkname = "pubmed_pubmed_citedin"),
    max_retries = max_retries + 2, delay = delay * 15,
    label = paste("entrez_link for", PubMedID)
  )
  if(is.null(src)){
    print(paste("Could not retrieve link data for", PubMedID, "- skipping", sep = " "))
    return(NULL)
  }

  print(paste("Extract list of publications citing", PubMedID, sep = " "))
  #Extract PubMedIDs of pubs citing our target Pub (= vect of pubmedIDs)
  citedIn <- src$links$pubmed_pubmed_citedin
  print(paste(PubMedID, "is cited in", length(citedIn), "publications", sep = " "))

  # If nothing cites this PubMedID, there is nothing left to fetch.
  # (Looping over 1:length(citedIn) when length==0 would run for i = 1 and 0
  # and crash on citedIn[1]/citedIn[0], so we bail out early instead.)
  if(is.null(citedIn) || length(citedIn) == 0){
    print(paste("No citing publications found for", PubMedID, "- skipping", sep = " "))
    return(NULL)
  }

  # Fetch citing records in batches instead of one request per PMID.
  # Firing one entrez_fetch() call per citation means thousands of individual
  # HTTP requests for a well-cited paper - NCBI throttles/blocks that pattern,
  # and it also means many more chances for a dropped connection. Grabbing
  # ~50 records per request cuts a 2000-citation run down to ~40 requests.
  chunks <- split(citedIn, ceiling(seq_along(citedIn) / batch_size))

  OUT <- NULL
  print(paste("Download data on publications citing", PubMedID, "in", length(chunks),
              "batch(es) of up to", batch_size, "records", sep = " "))
  pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)

  for(b in seq_along(chunks)){
    ids <- chunks[[b]]

    # Fetch one batch of XML records and parse it directly to a data frame in
    # one step. table_articles_byAuth() is the current exported, documented
    # easyPubMed function for this - the lower-level article_to_df()/
    # articles_to_list() helpers used in older examples are not exported in
    # this package version, so we go straight through the public API instead.
    # included_authors = "first" keeps exactly one row per article.
    batch_df <- .retry(
      function(){
        RefData <- rentrez::entrez_fetch(db = "pubmed", id = ids, rettype = "xml")
        easyPubMed::table_articles_byAuth(
          pubmed_data = RefData,
          included_authors = "first",
          max_chars = 500,
          autofill = FALSE,
          getKeywords = FALSE
        )
      },
      max_retries = max_retries + 2, delay = delay * 15,
      label = paste("batch", b, "of", length(chunks), "for", PubMedID)
    )

    if(!is.null(batch_df)) OUT <- rbind(OUT, batch_df)

    # Be polite to NCBI's E-utilities rate limit (3 req/sec without an API key,
    # 10 req/sec with one set via rentrez::set_entrez_key()).
    Sys.sleep(delay)

    # update progress bar
    setTxtProgressBar(pb, b)
  }
  close(pb)

  if(is.null(OUT) || nrow(OUT) == 0){
    print(paste("No publication data could be retrieved for citations of", PubMedID, sep = " "))
    return(NULL)
  }

  #Add col with TargetPubMed ID
  TargetPubMed <- rep(PubMedID, nrow(OUT))

  #FINAL dataset
  FINAL <- cbind(TargetPubMed, OUT)

  FINAL
}
