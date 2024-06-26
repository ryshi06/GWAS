
## should replace this with a method!
.getMetadata <- function(object, varname) {
  varMetadata(object@snpAnnot)[varname, "labelDescription"]
}

.setFilter <- function(genoData, filter.cols) {
  filt.df <- data.frame(getSnpVariable(genoData, filter.cols))
  names(filt.df) <- filter.cols
  filter <- rep("", nsnp(genoData))
  meta <- character()
  for (c in filter.cols) {
    filter[filt.df[[c]]] <- paste(filter[filt.df[[c]]], c, sep=";")
    meta[c] <- paste0('##FILTER=<ID=', c,
                      ',Description="', .getMetadata(genoData, c), '">')
  }
  filter <- sub("^;", "", filter)
  filter[filter == ""] <- "PASS"
  list(filter=filter, meta=meta)
}

.setInfo <- function(genoData, info.cols) {
  info.df <- data.frame(getSnpVariable(genoData, info.cols))
  names(info.df) <- info.cols
  info <- rep("", nsnp(genoData))
  meta <- character()
  for (c in info.cols) {
    if (is.logical(info.df[[c]])) {
      info[info.df[[c]]] <- paste(info[info.df[[c]]], c, sep=";")
    } else {
      info <- paste(info, paste0(c, "=", info.df[[c]]), sep=";")
    }
    type.map <- c(integer="Integer", numeric="Float", character="String",
                  factor="String", logical="Flag")
    meta[c] <- paste0('##INFO=<ID=', c,
                      ',Number=', ifelse(is.logical(info.df[[c]]), 0, 1),
                      ',Type=', type.map[class(info.df[[c]])],
                      ',Description="', .getMetadata(genoData, c), '">')
  }
  info <- sub("^;", "", info)
  info[info == ""] <- "."
  list(info=info, meta=meta)
}

vcfWrite <- function(genoData, vcf.file="out.vcf", sample.col="scanID",
                     id.col="snpID", qual.col=NULL, filter.cols=NULL,
                     info.cols=NULL, scan.exclude=NULL, snp.exclude=NULL,
                     scan.order=NULL,
                     ref.allele=NULL, block.size=1000, verbose=TRUE) {
  ## fixed fields
  chrom <- getChromosome(genoData, char=TRUE)
  pos <- getPosition(genoData)
  id <- getSnpVariable(genoData, id.col)
  ## check for missing values in id
  id <- as.character(id)
  id[is.na(id) | id == ""] <- "."
  if (!is.null(ref.allele)) {
    # stopifnot(length(ref.allele) == nsnp(genoData))
    # stopifnot(all(ref.allele %in% c("A", "B")))
    # a <- getAlleleA(genoData)
    # b <- getAlleleB(genoData)
    ref <- ref.allele
    alt <- alt.allele
  } else {
    ref <- "."
    alt <- "."
  }
  if (is.null(qual.col)) {
    qual <- rep(".", nsnp(genoData))
  } else {
    qual <- getSnpVariable(genoData, qual.col)
  }
  if (is.null(filter.cols)) {
    filter <- rep(".", nsnp(genoData))
    filt.meta <- character()
  } else {
    filt.list <- .setFilter(genoData, filter.cols)
    filter <- filt.list[["filter"]]
    filt.meta <- filt.list[["meta"]]
    rm(filt.list)
  }
  if (is.null(info.cols)) {
    info <- rep(".", nsnp(genoData))
    info.meta <- character()
  } else {
    info.list <- .setInfo(genoData, info.cols)
    info <- info.list[["info"]]
    info.meta <- info.list[["meta"]]
    rm(info.list)
  }
  format <- rep("GT", nsnp(genoData))
  fixed <- cbind(chrom, pos=as.character(pos), id=id,
                 ref, alt, qual=as.character(qual), filter,
                 info, format)
  ## subset with snp.exclude
  if (!is.null(snp.exclude)) {
    snp.index <- !(getSnpID(genoData) %in% snp.exclude)
  } else {
    snp.index <- rep(TRUE, nsnp(genoData))
  }
  
  ## open output file
  con <- file(vcf.file, "w")
  
  ## write metadata
  meta <- c('##fileformat=VCFv4.1',
            paste0('##fileDate=', Sys.Date()),
            '##source=GWASTools::vcfWrite()',
            filt.meta,
            info.meta,
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  writeLines(meta, con)
  
  ## get samples
  samples <- getScanVariable(genoData, sample.col)
  scanID <- getScanID(genoData)
  if (is.null(scan.order)) {
    scan.order <- scanID
  } else {
    stopifnot(all(scan.order %in% scanID))
  }
  ## subset with scan.exclude
  if (!is.null(scan.exclude)) {
    scan.order <- setdiff(scan.order, scan.exclude)
  }
  ## put samples in new order if requested
  scan.index <- match(scan.order, scanID)
  
  ## write header
  hdr <- paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                 samples[scan.index]), collapse="\t")
  writeLines(hdr, con)
  
  ## fwrite cannot write to an open connection
  close(con)
  
  ## write genotypes in blocks
  nblocks <- ceiling(nsnp(genoData) / block.size)
  for (i in 1:nblocks) {
    start <- (i-1)*block.size + 1
    end <- min(i*(block.size), nsnp(genoData))
    count <- end - start + 1
    n <- sum(snp.index[start:end])
    if (verbose) message("Block ", i, " of ", nblocks, "... ", n, " SNPs")
    
    if (n > 0) {
      geno <- getGenotype(genoData, snp=c(start,count), scan=c(1,-1))
      ## switch allele coding if ref.allele is B
      # geno[ref.allele[start:end] == "B",] <- 2 - geno[ref.allele[start:end] == "B",]
      geno <- geno[snp.index[start:end], scan.index, drop=FALSE]
      geno[is.na(geno)] <- "./."
      geno[geno == 2] <- "0/0"
      geno[geno == 1] <- "0/1"
      geno[geno == 0] <- "1/1"
      
      out <- cbind(fixed[(start:end)[snp.index[start:end]],,drop=FALSE], geno)
      fwrite(out, vcf.file, append=TRUE, quote=FALSE, sep="\t",
             row.names=FALSE, col.names=FALSE)
    }
  }
}