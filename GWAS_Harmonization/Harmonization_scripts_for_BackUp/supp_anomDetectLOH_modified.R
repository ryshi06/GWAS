library(DNAcopy)

getbdry <- function(eta, nperm, max.ones, tol= 1e-2) {
  bdry <- rep(0, max.ones*(max.ones+1)/2)
  zz <- .Fortran("getbdry",
                 as.double(eta),
                 as.integer(max.ones),
                 as.integer(nperm),
                 as.integer(max.ones*(max.ones+1)/2),
                 bdry=as.integer(bdry),
                 etastr=double(max.ones),
                 as.double(tol),
                 PACKAGE="DNAcopy")
  #  list("eta.star"=zz$etastr, "boundary"=zz$bdry)
  zz$bdry
}

### LOH raw determination ###
########## SUB-functions ##########

### Find RUNS ##########
.LOHruns<-function(index,geno,eleft,eright,run.size,inter.size) {
  ## index are eligible snp indices for a sample/chrom
  ## geno is genotypes corresponding to eligible snp's
  ## eleft and eright:endpoints of interval of interest given as snp indices
  ## run.size is minimum number of wanted markers to declare a run
  ## inter.size is number of unwanted markers allowed to interrupt a run
  
  if(eright<eleft) stop("non-existent interval given to LOHruns")
  indices<-index>=eleft & index<=eright
  if(sum(indices)<run.size){out<-NULL; return(out)}
  
  index1<-index[indices]
  geno1<-geno[indices]
  
  w<-rle(as.vector(geno1))
  #w<-rle(  ) has w[[1]] = lengths of runs, w[[2]] corresponding values
  vals<-w[[2]];lngs<-w[[1]]
  r0<-vals
  
  rlen<-length(r0)
  if(rlen==1){ if(r0==0){left<-index1[1];right<-index1[length(index1)]
  out<-data.frame(left,right)
  names(out)<-c("left","right")
  return(out)}  else {out<-NULL;return(out)}  }
  
  ##establish initial positions of alternating runs
  in.pos<-c(1,rep(NA,rlen-1))
  for(i in 2:rlen) in.pos[i]<-in.pos[i-1]+lngs[i-1]
  
  ##merging homoz intervals if separated by inter.size no. of hets
  # assuming sum of lengths of homo intervals on either side meets run.size criterion
  
  tpos<-which(r0==1&lngs<=inter.size) #identify small runs of hets
  smf<-which(r0==0 & lngs<run.size) #identify 'small' runs of homos
  if(length(tpos)!=0){
    if(tpos[1]==1) {
      if(lngs[2]>=run.size){r0[1]<-0}
      tpos<-tpos[-1]}
    if(length(tpos)!=0){
      if(tpos[length(tpos)]==rlen) {tpos<-tpos[-length(tpos)]; if(lngs[rlen-1]>=run.size) {r0[rlen]<-0}}
      if(length(tpos)!=0){
        for(k in tpos){if((lngs[k-1]+lngs[k+1])>=run.size) {r0[k]<-0}}
      }}  }
  
  ##want smaller runs of 0s to become runs of 'hets' but not if they
  #are part of a combined run
  run2<-rle(r0)
  vals2<-run2[[2]];lngs2<-run2[[1]]
  w0<-which(vals2==0 & lngs2>1)
  if(length(w0)==0){r0[smf]<-1} else {
    index.set<-NULL
    for(j in 1:length(w0)){ if(w0[j]==1){start<-1} else {start<-sum(lngs2[1:(w0[j]-1)])+1}
      end<-sum(lngs2[1:w0[j]])
      index.set<-c(index.set,start:end)}
    smf.use<-setdiff(smf,intersect(smf,index.set))
    r0[smf.use]<-1}
  
  ##after merging some of the initial runs, we get modified listing of r0 vals
  ## look for runs here; e.g. if run of two 0s that means that initially
  #we had a run of homo with small run of hets - rle length of 2 then indicates
  #putting those two original runs together as one run of homo 
  new.rle<-rle(r0)
  nvals<-new.rle[[2]]
  nlens<-new.rle[[1]]  ##indicates how many of original runs to put together
  if(length(nvals)==1){
    if(nvals!=0){out<-NULL;return(out)} else {
      left<-index1[1];right<-index1[length(index1)]
      out<-data.frame(left,right)
      names(out)<-c("left","right")
      return(out)} }
  
  newt<-which(nvals==0) #newt could be empty if originally there were no long homo runs
  if(length(newt)==0){out<-NULL;return(out)}
  left<-NULL
  right<-NULL
  ### if newt indicates runs of 0s, change initial and end positions of 
  #runs of homozy accordingly
  if(newt[1]==1){left<-c(left,1);right<-c(right,in.pos[nlens[1]+1]-1);newt<-newt[-1]}
  for(k in newt){ind<-sum(nlens[1:(k-1)]);left<-c(left,in.pos[ind+1])
  kl<-length(newt);kk<-newt[kl]
  if((ind+1+nlens[kk])<=length(in.pos)){
    right<-c(right,in.pos[ind+1+nlens[k]]-1)} else {right<-c(right,length(index1))}}
  ##right and left positions are indices of geno1, indices of index1
  if(length(right)==0|length(left)==0){out<-NULL;return(out)}
  out<-data.frame(index1[left],index1[right])
  names(out)<-c("left","right")
  return(out)
}#end function
################

########## Find BASE info #############
.LOHbase<-function(anoms,index,lrr,min.lrr.num) {
  ## anoms are anomalies found (BAF as well as potential LOH)  
  # for a given sample/chrom with left, right given as snp indices
  ## index,lrr,pos are for a given sample/chrom
  #min.lrr.num - minimum number of lrr pts to be considered/processed
  
  if(dim(anoms)[1]==0){newright<-index[length(index)]
  newleft<-index[1]} else {anoms<-anoms[order(anoms$left),]
  newright<-c(anoms$left,index[length(index)])
  newleft<-c(index[1],anoms$right)
  }
  
  lrr.pts<-NULL
  for(k in 1:length(newright)){ 
    
    int.ind<-index>=newleft[k]& index<=newright[k] 
    if(sum(int.ind)<min.lrr.num) next
    
    y<-lrr[int.ind]
    lrr.pts<-c(lrr.pts,y)
  }
  if(length(lrr.pts)<min.lrr.num) {out<-data.frame(NA,NA,NA,NA)} else{
    
    chmad<-mad(lrr.pts,na.rm=TRUE)
    chmed<-median(lrr.pts,na.rm=TRUE)
    chmn<-mean(lrr.pts,na.rm=TRUE)
    chsd<-sd(lrr.pts,na.rm=TRUE)
    out<-data.frame(chmad,chmed,chmn,chsd)
  }
  names(out)<-c("chrom.nonanom.mad","chrom.nonanom.median","chrom.nonanom.mean","chrom.nonanom.sd")
  return(out)  }
################################## 
########### MAIN FUNCTION ##############

.LOHfind<-function(snum,ch,geno,index,lrr,chr,segs,
                   ansch,run.size=50,inter.size=4,min.lrr.num=20 ) {
  
  #snum - current sample number
  #ch - current chromosome 
  #geno - genotypes for selected snps for that chrom
  #index - snp indices for selected snps for that chrom
  #lrr - LogRRatio values for selected snps for that chrom
  #chr - vector corresponding to given ch (same length as index, lrr)
  #segs - segments for given snum/ch from DNAcopy
  #ansch - previously found anomalies (e.g. BAF) to remove before finding runs
  #run.size - number of markers needed to declare a run
  #inter.size - number of consecutive hets allowed in a run
  #min.lrr.num - minimum number of lrr pts to be considered
  
  if(!is.element(class(ansch),"data.frame")) stop(paste("known anom input for sample ",snum," chromosome ",ch," to LOHfind needs to be data.frame",sep=""))
  if(!is.element(class(segs),"data.frame")) stop(paste("segs input to LOHfind for sample ",snum," chromosome ",ch," needs to be data.frame",sep=""))
  if(any(ansch$left>=ansch$right))stop(paste("some left >= right for known anoms for sample ",snum," chromosome ",ch,sep=""))
  ###
  num.segs<-dim(segs)[1]
  annum<-dim(ansch)[1]
  
  ####### find RUNS ###########
  ## find homozygous runs in intervals where there are no known anomalies
  ## (i.e. take out intervals determined by known anomalies)
  RUNS<-NULL
  if(annum!=0){
    LF<-c(index[1],ansch$right)
    RT<-c(ansch$left,index[length(index)])} else {
      LF<-index[1];RT<-index[length(index)]    }
  
  N<-length(LF)
  for(j in 1:N){ eleft<-LF[j];eright<-RT[j]
  outrun<-.LOHruns(index,geno,eleft,eright,run.size,inter.size)  #looking for runs in regions outside anoms
  RUNS<-rbind(RUNS,outrun)  }  #RUNS$right and $left are snp indices
  #################### find base info ########
  if(is.null(RUNS)) { 
    base<-.LOHbase(ansch,index,lrr,min.lrr.num)
    base$num.runs<-0
    base$num.segs<-num.segs
    base$scanID<-snum
    base$chrom<-ch
    
    if(!is.na(base$chrom.nonanom.mean) & !is.na(base$chrom.nonanom.sd)) {
      segs$sd.fac<-(segs$seg.mean-base$chrom.nonanom.mean)/base$chrom.nonanom.sd } else {segs$sd.fac<-NA}
    rout<-list(RUNS,base,segs)
    names(rout)<-c("RUNS","base.info","segments")
    return(rout)  }
  ###########
  num.runs<-dim(RUNS)[1]
  RUNS$scanID<-snum
  RUNS$chrom<-ch
  
  ## Base line info: chrom.mad, etc are computed for 'non.anom' region####
  ##  'non.anom' region excludes BAF anomalies and homoz runs (since these are potential anoms)####
  anoms1<-ansch[,c("left","right")]
  anoms2<-RUNS[,c("left","right")]
  anoms<-rbind(anoms1,anoms2)
  anoms<-anoms[order(anoms$left),]
  base<-.LOHbase(anoms,index,lrr,min.lrr.num)
  base$num.runs<-num.runs
  base$num.segs<-num.segs
  base$scanID<-snum
  base$chrom<-ch
  if(!is.na(base$chrom.nonanom.mean) & !is.na(base$chrom.nonanom.sd)) {
    segs$sd.fac<-(segs$seg.mean-base$chrom.nonanom.mean)/base$chrom.nonanom.sd } else {segs$sd.fac<-NA}
  rout<-list(RUNS,base,segs)
  names(rout)<-c("RUNS","base.info","segments")
  return(rout)  
  
} #end of function

############ Find anomalous segments from homozygous runs ##############
#### subfunctions to merge consecutive intervals passing filter #########
.mmerge<-function(w,segs){
  N<-length(w)
  if(N<=1) {flag<-0;tmp2<-segs[w,c("left","right")]
  out<-list(flag,tmp2);
  names(out)<-c("flag","anoms"); return(out)}
  wdiff<-sapply(2:N,function(w,i){w[i]-w[i-1]},w=w)
  rn<-rle(wdiff)
  vals<-rn[[2]]
  lng<-rn[[1]]
  wone<-which(vals==1)
  if(length(wone)==0) { flag<-0;tmp2<-segs[w,c("left","right")]; out<-list(flag,tmp2);
  names(out)<-c("flag","anoms"); return(out)}
  
  rlen<-length(vals)
  if(rlen<2){init.pos<-1} else{
    init.pos<-c(1,rep(NA,rlen-1))
    for(i in 2:rlen) init.pos[i]<-init.pos[i-1]+lng[i-1]  }
  
  n1<-length(wone)
  comb.list<-list()
  for(j in 1:n1){
    comb.list[[j]]<-init.pos[wone[j]]:(init.pos[wone[j]]+lng[wone[j]]) 
    ## comb.list are positions of w
    ##we now need positions in list of segments
    comb.list[[j]]<-w[comb.list[[j]]]}
  
  mod.ind<-unlist(comb.list)
  an.comb<-NULL
  for(j in 1:n1) {
    ind<-comb.list[[j]]
    left<-min(segs$left[ind])
    right<-max(segs$right[ind])
    tmp<-data.frame(left,right)
    names(tmp)<-c("left","right")
    an.comb<-rbinddt(an.comb,tmp)   }
  ws<-setdiff(w,mod.ind)
  an.nomod<-segs[ws,c("left","right")]
  tmp2<-rbinddt(an.nomod,an.comb)
  flag<-1;out<-list(flag,tmp2)
  names(out)<-c("flag","anoms")
  return(out)   } #end function mmerge
################### main merge ##############
.runsMerge<-function(segs,sig) {
  ## output is new list of runs with merged endpoints ##
  
  wup<-which(segs$sd.fac>sig)
  wdown<-which(segs$sd.fac < (-sig) )
  del<-union(wup,wdown)
  if(length(del)==0){rest<-segs[,c("left","right")]} else {
    
    rest<-segs[-del,c("left","right")] }
  
  ### up ###
  tmpup<-NULL
  if(length(wup)!=0){
    resup<-.mmerge(wup,segs)
    tmpup<-resup$anoms  } # end if(wup!=0)
  ## down ## 
  tmpdown<-NULL
  if(length(wdown)!=0){
    resdown<-.mmerge(wdown,segs)
    tmpdown<-resdown$anoms
  } # end if(wdown!=0)
  #####
  flag<-FALSE
  if(length(wup)!=0) { if(resup$flag==1) flag<-TRUE}
  if(length(wdown)!=0) { if (resdown$flag==1) flag<-TRUE}
  out<-rbinddt(rest,tmpup,tmpdown)
  out<-out[order(out$left),]
  rout<-list(out,flag)
  names(rout)<-c("newsegs","flag")
  return(rout)
} #end runsMerge 
############################################################
############local mad #########################
.LOHlocalMad<-function(select,index,nonanom.index,lrr,length.factor,min.lrr.num){
  if(any(select$left>select$right))stop("some left endpts > right endpts in runs input dataframe")
  tm5<-NULL
  for(j in 1:dim(select)[1]){
    lt<-select$left[j]; rt<-select$right[j]
    int<-index>=select$left[j] & index<=select$right[j] 
    xx<-lrr[int]
    xxm<-median(xx,na.rm=TRUE)
    
    zf<-length.factor*sum(int)
    left.bdy<-index[1]
    right.bdy<-index[length(index)]
    nleft<-lt-zf; nright<-rt+zf
    nleft<-max(nleft,left.bdy);nright<-min(nright,right.bdy)
    int.test<-index>=nleft & index<=nright
    nonanom<-is.element(index,nonanom.index)
    int5<-int.test&nonanom
    yy<-lrr[int5]
    if(length(yy)< min.lrr.num){tm5<-c(tm5,NA)} else { 
      yym<-median(yy,na.rm=TRUE);yymad<-mad(yy,na.rm=TRUE)
      tm<-(xxm-yym)/yymad
      tm5<-c(tm5,tm)  }
    
  } #end of loop on j
  out<-select
  out$local<-tm5
  return(out)
} #end function

############# MAIN ########################
## find breakpoints in a homozygous run ##
##############

.LOHselectAnoms<-function(snum,ch,segs,RUNS,base,index,nonanom.index,lrr,
                          homodel.min.num=10,homodel.thresh=10,small.num=20,small.thresh=2.25, 
                          medium.num=50,medium.thresh=2,long.num=100,long.thresh=1.5, small.na.thresh=2.5,
                          length.factor=5,merge.fac=.85,min.lrr.num=20) {
  
  #segs: DNAcopy segments for current sample/chromsome
  #
  if(!all(is.element(class(RUNS),"data.frame"))) stop("RUNS needs to be a data.frame")
  if(!all(is.element(c("scanID","chrom","left","right"),names(segs)))) stop("some names of RUNS missing")
  if(!all(is.element(class(segs),"data.frame"))) stop("segs needs to be a data.frame")
  if(!all(is.element(c("scanID","chrom","left","right","sd.fac"),names(segs)))) stop("some names of segs missing")
  if(!all(is.element(class(base),"data.frame"))) stop("base needs to be a data.frame")
  if(!all(is.element(c("scanID","chrom","chrom.nonanom.median","chrom.nonanom.mad","chrom.nonanom.mean","chrom.nonanom.sd"),names(base)))) stop("some names of base missing")
  
  if(!all(RUNS$scanID==snum & RUNS$chrom==ch)) stop("RUNS not from same sample/chrom")
  if(!all(segs$scanID==snum & segs$chrom==ch)) stop("segs not from same sample/chrom")
  if(!all(base$scanID==snum & base$chrom==ch)) stop("base not from same sample/chrom")
  
  ##### merge segs if appropriate #####
  new<- .runsMerge(segs,merge.fac)
  FLAG<-new$flag
  select<-new$newsegs
  
  #### find overlap of segs with RUNS ##########
  if(dim(RUNS)[1]<1) stop("No runs to test")
  if(dim(segs)[1]<1)stop("No segment info")
  runsegs<-NULL
  for(i in 1:dim(RUNS)[1]) {
    for(k in 1:dim(select)[1] ){
      mxL<-max(c(RUNS$left[i],select$left[k]))
      mnR<-min(c(RUNS$right[i],select$right[k]))
      over<-mnR-mxL
      
      if(over<=0)next  #no overlap
      tp<-data.frame(mxL,mnR)
      names(tp)<-c("left","right")
      runsegs<-rbinddt(runsegs,tp)    
    } 
  }
  
  ## find stats for the new runs
  sd.fac<-NULL;mad.fac<-NULL;seg.median<-NULL
  seg.median<-NULL;seg.mean<-NULL;nm<-NULL
  for(k in 1:dim(runsegs)[1]){
    subind<-index>=runsegs$left[k] & index<=runsegs$right[k]
    nm<-c(nm,sum(subind))
    sublrr<-lrr[subind]
    seg.med<-median(sublrr,na.rm=TRUE)
    seg.median<-c(seg.median,seg.med)
    sgm<-mean(sublrr,na.rm=TRUE)
    seg.mean<-c(seg.mean,sgm)
    sdf<-(sgm-base$chrom.nonanom.mean)/base$chrom.nonanom.sd
    mdf<-(seg.med-base$chrom.nonanom.median)/base$chrom.nonanom.mad
    sd.fac<-c(sd.fac,sdf); mad.fac<-c(mad.fac,mdf) 
  }#end of k loop
  runsegs$num.mark<-nm
  runsegs$seg.median<-seg.median
  runsegs$seg.mean<-seg.mean
  runsegs$mad.fac<-mad.fac 
  runsegs$sd.fac<-sd.fac
  runsegs$scanID<-snum
  runsegs$chrom<-ch
  runsegs$chrom.nonanom.median<-base$chrom.nonanom.median
  runsegs$chrom.nonanom.mad<-base$chrom.nonanom.mad
  runsegs$chrom.nonanom.mean<-base$chrom.nonanom.mean
  runsegs$chrom.nonanom.sd<-base$chrom.nonanom.sd
  ### select ##
  
  num.segs<-dim(segs)[1]
  
  ## LOHlocalMad computes and adds the variable "local"
  select2<-.LOHlocalMad(runsegs,index,nonanom.index,lrr,length.factor,min.lrr.num)
  select2$num.segs<-num.segs
  
  lenrest<-dim(select2)[1]
  #will have exited before if there is nothing in select2
  
  ## apply filters
  jind<-NULL
  for(j in 1:lenrest){ 
    
    if(is.na(select2$local[j])) {md<-abs(select2$mad.fac[j]);minv<-md} else {
      md<-mean(c(abs(select2$mad.fac[j]),abs(select2$local[j])))
      minv<-min(c(abs(select2$mad.fac[j]),abs(select2$local[j])))  
    }
    nm<-select2$num.mark[j]
    if(nm>=homodel.min.num & abs(select2$mad.fac[j])>homodel.thresh){ jind<-c(jind,j);next}
    if(!is.na(select2$local[j])& nm>small.num & nm<=medium.num &md>small.thresh & minv>=medium.thresh){ jind<-c(jind,j);next}
    if(is.na(select2$local[j])& nm>small.num & nm<=medium.num & md>small.na.thresh) {jind<-c(jind,j);next}
    if(nm>medium.num & nm<=long.num &md>medium.thresh & minv>=long.thresh) {jind<-c(jind,j);next}
    if(nm>long.num & md>long.thresh & signif(minv,2)>=long.thresh) {jind<-c(jind,j);next}
    
  } #end j loop
  if(length(jind)!=0) {selt<-select2[jind,];selt<-selt[order(selt$left),]} else selt<-NULL
  outt<-list(select2,selt,FLAG)
  names(outt)<-c("raw.adjusted","filtered","merge.flag")
  return(outt)
} #end of function

anomDetectLOH<-function(intenData, genoData, scan.ids, chrom.ids, snp.ids,
                        known.anoms,
                        smooth=50,min.width=5,nperm=10000,alpha=.001,
                        run.size=50,inter.size=4,
                        homodel.min.num=10,homodel.thresh=10,
                        small.num=20,small.thresh=2.25, medium.num=50,medium.thresh=2,
                        long.num=100,long.thresh=1.5, small.na.thresh=2.5,
                        length.factor=5,merge.fac=.85,min.lrr.num=20,verbose=TRUE){
  
  ##scan.ids: samples to consider
  ##chrom.ids: chromosomes to consider, usually 1:23
  ##snp.ids: intid's for eligible snps
  ##known.anoms: data.frame of known anoms (usually from BAF); 
  #   must have "scanID","chrom","left","right" where left and right are snp indices
  
  ##run.size: number of markers to declare a 'homozygous' run (here homozygous includes missing) 
  ##inter.size: number of consecutive heterozygous markers allowed to interrupt a run
  ### detection of anomalies based on a chromsome-wide and local mad.fac thresholds
  #    where mad.fac is (segment median-nonanomalous median)/nonanom mad
  ##homodel.min.num: minimum number of markers to detect extreme diff in lrr (for homo deletion)
  ##homodel.thresh: threshold for 'mad.fac' to detect extreme diff in lrr
  ##small.num: minimum number of markers to detect anomaly (other than extreme)
  ##small.thresh:threshold for  number of markers between small.num and medium.num
  ##medium.num: number of markers threshold for 'medium' 
  ##medium.thresh: threshold for 'mad.fac' when number of markers is between medium num and long num
  ##long.num: number of markers threshold for 'long'
  ##long.thresh:threshold for 'mad.fac' when number of markers is bigger than long.num
  ##small.na.thresh: 
  #   chrom mad.fac threshold when between small.num and medium.num and no local mad.fac
  ##length.factor: 
  #   local mad.fac based on interval that is length.factor*(no. of markers in segment) on either side
  ##merge.fac: threshold used to merge original segmentation segments 
  ##min.lrr.num: if any 'nonanomalous' interval less than min.lrr.num,
  #  ignore this piece in finding overall nonanomalous unless is only piece left
  
  #### checks ####
  # check that intenData has LRR
  if (!hasLogRRatio(intenData)) stop("LogRRatio not found in intenData")
  
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
  
  # check that sex is present in annotation
  if (hasSex(intenData)) {
    sex <- getSex(intenData)
  } else if (hasSex(genoData)) {
    sex <- getSex(genoData)
  } else stop("sex not found in intenData or genoData")
  
  intid <- intenSnpID
  if (!all(is.element(snp.ids,intid))) stop("snp.ids has values not present in intenData")
  
  chrom <- getChromosome(intenData)
  if (!all(is.element(chrom.ids, chrom))) stop("chrom.ids has values not present in intenData (all values in chrom.ids should be integers)")
  
  sid <- intenScanID
  male <- sid[!is.na(sex) & sex == "M"]
  
  if(!is.element(class(known.anoms),"data.frame") | !all(is.element(c("scanID","chromosome","left.index","right.index"),names(known.anoms)))){
    stop("known.anoms input needs to be data.frame with variables including scanID, chromosome, left.index, right.index") }
  ####
  
  # internal functions require these names, so convert from package standard
  names(known.anoms)[names(known.anoms) == "chromosome"] <- "chrom"
  names(known.anoms)[names(known.anoms) == "left.index"] <- "left"
  names(known.anoms)[names(known.anoms) == "right.index"] <- "right"
  
  LOH.raw<-NULL;LOH.base.info<-NULL; LOH.filtered<-NULL
  LOH.segments<-NULL; LOH.merge<-NULL;LOH.raw.adjusted<-NULL
  ## compute parameter needed for DNAcopy
  max.ones<-floor(nperm*alpha)+1
  sbdry<-getbdry(.05,nperm,max.ones)
  
  sel<-is.element(intid,snp.ids)
  orindex<-which(sel)  #indices of eligible snp's
  orchr<-chrom[sel]
  for(snum in scan.ids){
    sindex <- which(is.element(sid, snum))
    if(length(sindex)==0) stop(paste("Sample ",snum, " does not exist",sep=""))
    GENO <- getGenotype(genoData, snp=c(1,-1), scan=c(sindex,1))
    olrr<-getLogRRatio(intenData, snp=c(1,-1), scan=c(sindex,1))[sel]  
    ogeno<-GENO[sel]
    ws<-!is.na(olrr)
    olrr<-olrr[ws]
    oindex<-orindex[ws]
    ogeno<-ogeno[ws]
    chr<-orchr[ws]
    
    ##homoz and missing have value 0
    whomo<-is.element(ogeno,c(0,2)) | is.na(ogeno)
    ogeno[whomo]<-0
    
    ####### DNAcopy Segmentation for given sample ########
    temp.CNA<-DNAcopy::CNA(as.vector(olrr),chr,oindex,data.type="logratio",sampleid=snum)
    temp.smooth<-DNAcopy::smooth.CNA(temp.CNA,smooth.region=smooth,outlier.SD.scale=4)  #smooth.region,outlier.SD.scale=default
    temp.segment<-segment(temp.smooth,alpha=alpha,sbdry=sbdry,p.method="h",min.width=min.width,nperm=nperm,undo.splits="sdundo",undo.SD=1,verbose=as.integer(verbose))
    segments<-temp.segment$out
    if(dim(segments)[1]<1) next
    segments$ID<-rep(snum,length(segments$ID))
    #$loc.start and $loc.end are indices of snp
    names(segments)<-c("scanID","chrom","left","right","num.mark","seg.mean")
    
    
    for(ch in chrom.ids){
      if(ch==XchromCode(intenData) & is.element(snum,male)) next
      
      wc<-chr==ch 
      geno<-ogeno[wc]
      index<-oindex[wc]
      lrr<-olrr[wc]   
      chrr<-chr[wc]
      
      ansch<-known.anoms[known.anoms$scanID==snum & known.anoms$chrom==ch,]
      segs<-segments[segments$chrom==ch,]
      
      ###### find homozygous runs and base info for current sample/chrom
      out<-.LOHfind(snum,ch,geno,index,lrr,chrr,segs,ansch,
                    run.size,inter.size,min.lrr.num )
      
      base.snch<-out$base.info
      base.snch$sex<-sex[sindex]
      RUNS.snch<-out$RUNS
      segs.snch<-out$segments
      
      LOH.base.info<-rbinddt(LOH.base.info,base.snch)
      LOH.segments<-rbinddt(LOH.segments,segs.snch)
      LOH.raw<-rbinddt(LOH.raw,RUNS.snch)
      
      #### BEGIN filtering process  #############
      
      ####  Special Cases ###########  
      
      if(base.snch$num.runs==0) next  #no anomalies
      
      if(is.na(base.snch$chrom.nonanom.mad) & is.element(base.snch$num.runs,c(1,2)) & dim(ansch)[1]==0) {
        # happens if runs take up most of chromosome with no other known anoms, i.e. probably whole chrom LOH
        seg.median<-NULL;seg.mean<-NULL;nm<-NULL
        for(k in 1:dim(RUNS.snch)[1]){
          subind<-index>=RUNS.snch$left[k] & index<=RUNS.snch$right[k]         
          nm<-c(nm,sum(subind))
          sublrr<-lrr[subind]
          seg.med<-median(sublrr,na.rm=TRUE)
          seg.median<-c(seg.median,seg.med)
          seg.mean<-c(seg.mean,mean(sublrr,na.rm=TRUE))
        }#end of k loop
        RUNS.snch$num.mark<-nm
        RUNS.snch$seg.median<-seg.median
        RUNS.snch$seg.mean<-seg.mean
        RUNS.snch$sd.fac<-NA
        RUNS.snch$mad.fac<-NA
        RUNS.snch$local<-NA
        RUNS.snch$num.segs<-base.snch$num.segs
        RUNS.snch$chrom.nonanom.mad<-NA
        RUNS.snch$chrom.nonanom.median<-NA
        RUNS.snch$chrom.nonanom.mean<-NA
        RUNS.snch$chrom.nonanom.sd<-NA
        RUNS.snch$sex<-sex[sindex] 
        LOH.raw.adjusted<-rbinddt(LOH.raw.adjusted,RUNS.snch)
        LOH.filtered<-rbinddt(LOH.filtered,RUNS.snch)
        next
      }
      
      if(base.snch$num.segs==1) next
      #if no segmentation and previous situation not occur, there are no anomalies
      if(is.na(base.snch$chrom.nonanom.mad)) next  
      #no base left to compare with;already singled out whole chrom possibiltiy; too many pieces: no anoms
      
      ##### end Special cases###
      ### Get info needed for filtering process
      ## find indices for 'nonanom' snps: not in known anoms nor in any found runs (which are potential anoms) 
      baf.del<-NULL
      if(dim(ansch)[1]!=0){
        for(j in 1:dim(ansch)[1]) {
          int<-index[index>=ansch$left[j]& index<=ansch$right[j]]
          baf.del<-union(baf.del,int)
        }
      }  #index values
      
      runs.del<-NULL             #delete all found runs or anoms to determine "nonanom" base
      for(i in 1:dim(RUNS.snch)[1]){
        int<-index[index>=RUNS.snch$left[i] & index<=RUNS.snch$right[i]]
        runs.del<-union(runs.del,int) 
      }    
      
      possible.anom.index<-union(baf.del,runs.del)
      #by previous checks, will have at least some runs and runs+known anoms will not take up whole chrom
      nonanom.index<-setdiff(index,possible.anom.index)
      
      #### final FILTERING #########
      outt<-.LOHselectAnoms(snum,ch,segs.snch,RUNS.snch,base.snch,index,nonanom.index,lrr,
                            homodel.min.num,homodel.thresh,small.num,small.thresh, medium.num,medium.thresh,
                            long.num,long.thresh, small.na.thresh,length.factor,merge.fac,min.lrr.num)
      raw.adj<-outt$raw.adjusted
      if(!is.null(raw.adj)) raw.adj$sex<-sex[sindex]  
      filtered<-outt$filtered
      if(!is.null(filtered)) filtered$sex<-sex[sindex] 
      LOH.raw.adjusted<-rbinddt(LOH.raw.adjusted,raw.adj)
      LOH.filtered<-rbinddt(LOH.filtered, filtered)
      if(outt$merge.flag){
        tmp<-data.frame(snum,ch);names(tmp)<-c("scanID","chrom")
        LOH.merge<-rbinddt(LOH.merge,tmp)
      }
      
    } #end chrom loop
  } #end sample loop
  
  
  colm<-c("scanID","chrom","left","right","num.mark","seg.median","seg.mean",
          "mad.fac","sd.fac","local","num.segs","chrom.nonanom.mad","chrom.nonanom.median","chrom.nonanom.mean","chrom.nonanom.sd","sex")
  
  if(!is.null(LOH.raw.adjusted)){
    LOH.raw.adjusted<-LOH.raw.adjusted[,colm]
  }
  
  if(nrow(LOH.filtered) != 0){
    LOH.filtered<-LOH.filtered[order(LOH.filtered$scanID,LOH.filtered$chrom,LOH.filtered$left),]
    LOH.filtered<-LOH.filtered[,colm] 
  }
  
  outr<-list(LOH.raw,LOH.raw.adjusted,LOH.filtered,LOH.base.info,LOH.segments,LOH.merge)
  names(outr)<-c("raw","raw.adjusted","filtered","base.info","segments","merge")
  
  # convert back to package standard names
  for (i in 1:length(outr)) {
    if ("chrom" %in% names(outr[[i]]))
      names(outr[[i]])[names(outr[[i]]) == "chrom"] <- "chromosome"
    if ("left" %in% names(outr[[i]]))
      names(outr[[i]])[names(outr[[i]]) == "left"] <- "left.index"
    if ("right" %in% names(outr[[i]]))
      names(outr[[i]])[names(outr[[i]]) == "right"] <- "right.index"
  }
  
  # convert index to base position
  for (i in 1:length(outr)) {
    if (!is.null(outr[[i]]) & ("left.index" %in% names(outr[[i]]))) {
      outr[[i]]$left.base <- getPosition(intenData, index=outr[[i]]$left.index)
      outr[[i]]$right.base <- getPosition(intenData, index=outr[[i]]$right.index)
    }
  }
  
  return(outr)
}