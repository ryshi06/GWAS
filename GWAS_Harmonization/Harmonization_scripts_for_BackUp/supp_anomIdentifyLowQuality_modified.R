
#load an Rdata file and return the object
#(using load directly returns only the object name)
getobj <- function(Rdata) {
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

#save an R object obj as name in an Rdata file called path/name.RData
saveas <- function(obj, name, path=".") {
  assign(name, obj)
  if (grepl("[.]RData$", name, ignore.case=TRUE) |
      grepl("[.]rda$", name, ignore.case=TRUE)) {
    file <- file.path(path, name)
  }
  else {
    file <- file.path(path, paste(name, "RData", sep="."))
  }
  save(list=name, file=file)
}

#returns TRUE if x and y exist and all elements are equal
#if there are NA values, returns TRUE if is.na(x) == is.na(y) and
#all other elements are equal, otherwise returns FALSE
#retuns FALSE if x or y (but not both) is NULL
allequal <- function(x, y) {
  if ((is.null(x) & !is.null(y)) | (!is.null(x) & is.null(y))) return(FALSE)
  
  if (!all(is.na(x) == is.na(y))) return(FALSE)
  
  return(all(x[!is.na(x)] == y[!is.na(y)]))
}

#####
# Read and write the first n lines of a file
#####
readWriteFirst <- function (filein, fileout, n)
{
  # read first n lines of filein and write them to fileout, where filein and fileout are file names
  incon <- file(filein,"r")  # open connection for reading
  outcon <- file(fileout,"w")  # open for writing
  for(i in 1:n) 
    write(readLines(incon, n=1), file=outcon, append=TRUE)
  close(incon)
  close(outcon)
}

# faster rbind using data.table
rbinddt <- function(...) {
  as.data.frame(data.table::rbindlist(list(...)))
}


###
.sampCharac<-function(snp,med.sd,seg.info){
  #find segmentation factors and record along with sd info for BAF and LOH
  #segmentation factors are ratio of number of segments to number of eligible markers
  #  found for each chromosome 1-23 and one factor for all autosomes
  #med.sd is dataframe with "scanID" and "med.sd" (median sd across autosomes)
  #outputs data.frame with variables
  #   c("scanID","chrX.num.segs","chrX.fac","max.autosome","max.auto.fac","max.auto.num.segs","fac.all.auto","med.sd")
  
  if(!is.element(class(med.sd),"data.frame") | !all(is.element(c("scanID","med.sd"),names(med.sd)))) stop("incorrect class or names for med.sd")
  if(!is.element(class(seg.info),"data.frame") | !all(is.element(c("scanID","chromosome","num.segs"),names(seg.info)))) stop("incorrect class or names for seg.info") 
  if(!is.element(class(snp),"SnpAnnotationDataFrame"))  stop("snp class incorrect")
  if(!hasVariable(snp, "eligible")
     |!is.element(class(snp$eligible),"logical") ) stop("snp variables incorrect")
  ##################################################
  
  snp.chrom <- getChromosome(snp)
  
  num.elig<-sum(snp$eligible & is.element(snp.chrom, autosomeCode(snp)))
  ##needed when considering factor over autosomes
  #######################
  ##eligible marker counts for each chromosome
  ch.marker.info<-NULL
  for(i in unique(snp.chrom)) {w<-snp.chrom==i& snp$eligible
  m<-sum(w)
  tmp<-data.frame(i,m)
  names(tmp)<-c("chromosome","marker.length")
  ch.marker.info<-rbinddt(ch.marker.info,tmp)
  }
  #############
  del<-seg.info$chromosome!=XYchromCode(snp)
  seg.info<-seg.info[del,]
  if(dim(seg.info)[1]==0) stop("Error: no seg info")
  
  smp<-unique(seg.info$scanID) 
  sampchr<-NULL
  for(snum in smp){ 
    nums<-seg.info[seg.info$scanID==snum,] 
    bfsd<-med.sd$med.sd[med.sd$scanID==snum]
    facs<-NULL
    for(k in 1:dim(nums)[1]) {
      ch<-nums$chromosome[k]
      ml<-ch.marker.info$marker.length[ch.marker.info$chromosome==ch]
      fac<-nums$num.segs[k]/ml
      facs<-c(facs,fac)
    } 
    nums$facs<-facs
    wX<-which(nums$chromosome==XchromCode(snp))
    if(length(wX)!=0) {mfX<-nums$facs[nums$chromosome==XchromCode(snp)];nsX<-nums$num.segs[nums$chromosome==XchromCode(snp)]} else {
      mfX<-NA;nsX<-NA}
    nums2<-nums[nums$chromosome!=XchromCode(snp),]  
    mfauto<-max(nums2$facs)
    wauto<-which(nums2$facs==mfauto)
    mfach<-nums2$chromosome[wauto]
    anseg<-nums2$num.segs[wauto]
    fac.all<-sum(nums2$num.segs)/num.elig
    tch<-length(which(nums$num.segs>1))
    tmp<-data.frame(snum,nsX,mfX,mfach,mfauto,anseg,tch,fac.all,bfsd)
    names(tmp)<-c("scanID","chrX.num.segs","chrX.fac","max.autosome","max.auto.fac","max.auto.num.segs","num.ch.segd","fac.all.auto","med.sd")
    sampchr<-rbinddt(sampchr,tmp)
  }
  return(sampchr)
} #end function

## identify messy samples
.identifyLowQual<-function(samp.info,sd.thresh,sng.seg.thresh,auto.seg.thresh){
  req<-c("scanID","chrX.num.segs","chrX.fac","max.autosome","max.auto.fac","max.auto.num.segs","fac.all.auto","med.sd")
  if(!is.element(class(samp.info),"data.frame") | 
     !all(is.element(req,names(samp.info)))) stop("class or names of samp.info incorrect")
  
  chbf1<-samp.info$med.sd>=sd.thresh 
  chbf2a<-samp.info$max.auto.fac>=sng.seg.thresh 
  chbf2b<-!is.na(samp.info$chrX.fac) &  samp.info$chrX.fac>=sng.seg.thresh 
  chbf3<-samp.info$fac.all.auto>=auto.seg.thresh 
  
  chsng<-chbf2a & !chbf1 &!chbf3 &!chbf2b
  chsngX<-chbf2b & !chbf1 & !chbf3 & !chbf2a
  chbfsd<-chbf1 & !chbf3  
  chbfseg<-chbf3 & !chbf1 
  chbfboth<-chbf3& chbf1
  if(sum(chbfsd)!=0) {
    bad.sd<-samp.info[chbfsd,]
    bad.sd$type<-"sd"} else bad.sd<-NULL
  if(sum(chsng)!=0) {bad.sng<-samp.info[chsng,]
  bad.sng$type<-"sng.seg"} else bad.sng<-NULL
  if(sum(chsngX)!=0){bad.sngX<-samp.info[chsngX,]
  bad.sngX$type<-"sng.seg.X"} else bad.sngX<-NULL
  if(sum(chbfseg)!=0){bad.auto<-samp.info[chbfseg,]
  bad.auto$type<-"auto.seg"} else bad.auto<-NULL
  if(sum(chbfboth)!=0) {bad.both<-samp.info[chbfboth,]
  bad.both$type<-"both.sd.seg"} else bad.both<-NULL
  bad<-rbinddt(bad.sd,bad.sng,bad.sngX,bad.auto,bad.both)
  # modified here: although all dataframe are NULL, after rbindlist(), 
  ## it will generate a dataframe with 0 columns and 0 rows, which will cause !(is.null(bad)) returns TRUE, 
  ## and further cause the error saying "Error in order(bad$scanID) : argument 1 is not a vector" 
  # if(!is.null(bad)) bad<-bad[order(bad$scanID),]
  if(nrow(bad) != 0) bad<-bad[order(bad$scanID),]
  return(bad)
}# end function

anomIdentifyLowQuality <- function(snp.annot, med.sd, seg.info,
                                   sd.thresh,sng.seg.thresh,auto.seg.thresh) {
  samp.info <- .sampCharac(snp.annot, med.sd, seg.info)
  low.qual <- .identifyLowQual(samp.info, sd.thresh, sng.seg.thresh, auto.seg.thresh)
  return(low.qual)
}

