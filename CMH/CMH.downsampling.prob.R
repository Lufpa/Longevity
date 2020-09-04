args = commandArgs(trailingOnly=T)
print (R.Version()$version.string)
print ("arguments: ")
print (args)
print ("start pipeline")

#load tables with alt and ref alleles for all cages to be tested in one CMH
#these tables only contain SNPs that overlap all cages

file <- read.table(paste(args[7],"/",args[1],sep=""),h=T,check.names=F)
fileref <- read.table(paste(args[7],"/",args[2],sep=""),h=T,check.names=F)

print(c("#SNPS", dim(file)[1], "#flies", dim(file)[2]))

samplematrix <- function(altmatrix, refmatrix) {
  # Create matrix of probabilites of reference
  probref <- refmatrix / (altmatrix + refmatrix)
  # Boolean vector where TRUE means we sample a reference allele
  sample <- sapply(probref, function(x) return(runif(1) <= x))
  # Reshape vector to matrix
  #dim(sample) <- dim(altmatrix)
  # Apply colnames
  #colnames(sample) <- colnames(altmatrix)
  # Convert to matrix of 0 and 1 and transpose to match old function
  return(1*sample)
}

#create the output table with one read per individual  
countsmatrix <- data.frame(matrix(ncol = nrow(file), nrow=ncol(file[,-1])))

#number of permutations  
perm=args[3]
perm <- as.numeric(perm)
#create the output tables for each cage where the total counts per cage will be stored
if (args[6] == "2reps") {
  perm_ref_a1 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_a1 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_a2 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_a2 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_at01 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_at01 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_at02 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_at02 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  
  ind <- read.table(paste(args[7],"indpercage_",args[4],"_", args[5],".txt",sep=""), h=T, row.names=1)
  
  for (p in 1:perm){
    print(c("permutation number ", p))
    start= Sys.time()
    countsmatrix <- samplematrix(file[,-1], fileref[,-1])

      
    perm_ref_a1[,p] <- apply(countsmatrix[,ind[1,1]:ind[1,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
     perm_alt_a1[,p] <- abs(perm_ref_a1[,p] - (apply(countsmatrix[,ind[1,1]:ind[1,2]], 1, function(x) sum(!is.na(x))))) #alt counts
     perm_ref_a2[,p] <- apply(countsmatrix[,ind[2,1]:ind[2,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
     perm_alt_a2[,p] <- abs(perm_ref_a2[,p] - (apply(countsmatrix[,ind[2,1]:ind[2,2]], 1, function(x) sum(!is.na(x))))) #alt counts
     perm_ref_at01[,p] <- apply(countsmatrix[,ind[3,1]:ind[3,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_at01[,p] <- abs(perm_ref_at01[,p] - (apply(countsmatrix[,ind[3,1]:ind[3,2]], 1, function(x) sum(!is.na(x))))) #alt counts
     perm_ref_at02[,p] <- apply(countsmatrix[,ind[4,1]:ind[4,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
     perm_alt_at02[,p] <- abs(perm_ref_at02[,p] - (apply(countsmatrix[,ind[4,1]:ind[4,2]], 1, function(x) sum(!is.na(x))))) #alt counts

    end=Sys.time()
    print (end-start)
  }
  
  ####write tables with counts per cage per SNP per permutation to be used eventually to run 
  ####one CMH for each permutation
  permfiles <- c("perm_ref_a1","perm_alt_a1","perm_ref_a2","perm_alt_a2","perm_ref_at01","perm_alt_at01","perm_ref_at02","perm_alt_at02")
  for (f in 1:length(permfiles)){
    ff <- permfiles[f]
    ta <- get(ff)
    write.table(cbind("site"=file$site,ta), paste("Permutation_outdir/",args[4],".",args[5],".",args[6],".",ff,".",args[8],".txt",sep=""), quote=F, row.names = F, col.names = T, sep="\t")
  }
  print("permutated counts per table were writen to Permutation_outdir")
  
  ####create the table with the average of counts derived from the permutations
  mean_a1 <- cbind("mean_ref"=apply(perm_ref_a1, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_a1,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_a2 <- cbind("mean_ref"=apply(perm_ref_a2, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_a2,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_at01 <- cbind("mean_ref"=apply(perm_ref_at01, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_at01,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_at02<- cbind("mean_ref"=apply(perm_ref_at02, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_at02,1,function(x) round(mean(x),0))) #mean ref and alt counts
  
  print("finished: permutations and averaging of counts")
  
  ####making sync file for popoolation
  #### chr pos ref a:t:c:g:N:del
  syncfile <- cbind.data.frame("chr"=file$site, "pos" = file$site, "N", 
                               "a1" = apply(mean_a1, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "a2"=apply(mean_a2, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "at01"=apply(mean_at01, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "at02"=apply(mean_at02, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")))
  
  
} else {
  
  perm_ref_a1 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_a1 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_a2 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_a2 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_at01 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_at01 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_at02 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_at02 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_d <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_d <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_ref_dt0 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  perm_alt_dt0 <- as.data.frame(matrix(ncol=perm, nrow = nrow(file)))
  
  #loop through all SNPs and select one read per individual, then add up the alt and ref reads to get counts
  #per individual per cage
  ind <- read.table(paste(args[7],"/indpercage/","indpercage_",args[4],"_", args[5],".txt",sep=""), h=T, row.names=1)
  for (p in 1:perm){
    print(c("permutation number ",p))
    start= Sys.time()
    countsmatrix <- samplematrix(file[,-1], fileref[,-1])
   

    perm_ref_a1[,p] <- apply(countsmatrix[,ind[1,1]:ind[1,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_a1[,p] <- abs(perm_ref_a1[,p] - (apply(countsmatrix[,ind[1,1]:ind[1,2]], 1, function(x) sum(!is.na(x))))) #alt counts
    perm_ref_a2[,p] <- apply(countsmatrix[,ind[2,1]:ind[2,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_a2[,p] <- abs(perm_ref_a2[,p] - (apply(countsmatrix[,ind[2,1]:ind[2,2]], 1, function(x) sum(!is.na(x))))) #alt counts
    perm_ref_at01[,p] <- apply(countsmatrix[,ind[3,1]:ind[3,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_at01[,p] <- abs(perm_ref_at01[,p] - (apply(countsmatrix[,ind[3,1]:ind[3,2]], 1, function(x) sum(!is.na(x))))) #alt counts
    perm_ref_at02[,p] <- apply(countsmatrix[,ind[4,1]:ind[4,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_at02[,p] <- abs(perm_ref_at02[,p] - (apply(countsmatrix[,ind[4,1]:ind[4,2]], 1, function(x) sum(!is.na(x))))) #alt counts
    perm_ref_dt0[,p] <- apply(countsmatrix [,ind[5,1]:ind[5,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_dt0[,p] <- abs(perm_ref_dt0[,p] - (apply(countsmatrix[,ind[5,1]:ind[5,2]], 1, function(x) sum(!is.na(x))))) #alt counts
    perm_ref_d[,p] <- apply(countsmatrix[,ind[6,1]:ind[6,2]], 1, function(x) sum(x, na.rm=T)) #ref counts
    perm_alt_d[,p] <- abs(perm_ref_d[,p] - (apply(countsmatrix[,ind[6,1]:ind[6,2]], 1, function(x) sum(!is.na(x))))) #alt counts
    
    end= Sys.time()
    print (end-start)
  }
  
  ####write tables with counts per cage per SNP per permutation to be used eventually to run 
  ####one CMH for each permutation
  permfiles <- c("perm_ref_a1","perm_alt_a1","perm_ref_a2","perm_alt_a2","perm_ref_at01","perm_alt_at01","perm_ref_at02","perm_alt_at02","perm_ref_d","perm_alt_d","perm_ref_dt0","perm_alt_dt0")
  for (f in 1:length(permfiles)){
    ff <- permfiles[f]
    ta <- get(ff)
    write.table(cbind("site"=file$site,ta), paste("Permutation_outdir/",args[4],".",args[5],".",args[6],".",ff,".",args[8],".txt",sep=""), quote=F, row.names = F, col.names = T, sep="\t")
  }
  
  print("permutated counts per table were writen to Permutation_outdir")
  
  ####create the table with the average of counts derived from the permutations
  
  mean_a1 <- cbind("mean_ref"=apply(perm_ref_a1, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_a1,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_a2 <- cbind("mean_ref"=apply(perm_ref_a2, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_a2,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_at01 <- cbind("mean_ref"=apply(perm_ref_at01, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_at01,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_at02<- cbind("mean_ref"=apply(perm_ref_at02, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_at02,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_d <- cbind("mean_ref"=apply(perm_ref_d, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_d,1,function(x) round(mean(x),0))) #mean ref and alt counts
  mean_dt0 <- cbind("mean_ref"=apply(perm_ref_dt0, 1, function(x) round(mean(x),0)), "mean_alt"=apply(perm_alt_dt0,1,function(x) round(mean(x),0))) #mean ref and alt counts
  print("finished: permutations and averaging of counts")
  
  
  ####making sync file for popoolation
  #### chr pos ref a:t:c:g:N:del
  syncfile <- cbind.data.frame("chr"=file$site, "pos" = file$site, "N",
                               "a1" = apply(mean_a1, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "a2"=apply(mean_a2, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "at01"=apply(mean_at01, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "at02"=apply(mean_at02, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "d"=apply(mean_d, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")),
                               "dt0"=apply(mean_dt0, 1, function(x) paste(x[1],":",x[2],":", 0,":",0,":",0,":",0, sep="")))
}



#site is chr:pos, separate chr and pos in two different columns
#chr
if (substr(syncfile[1,1],1,1) == 4){
  syncfile[,1] <- substr(syncfile[,1], 1,1)
}else {
  syncfile[,1] <- substr(syncfile[,1], 1,2)
}
#pos
syncfile[,2] <-as.character(syncfile[,2])
if (syncfile[1,1] == 4){
  substr(syncfile[,2],1,2) <- "  "
}else {
  substr(syncfile[,2], 1,3) <- "   "
}
syncfile <- syncfile[order(as.numeric(syncfile$pos), decreasing = F), ]  # popoolation NEEDS sorted coordinates
print ("check that syncfile looks OK")
head(syncfile); dim(syncfile)

####write sync table to be used as input in Popoolation
### This table's order is a1, a2, t01, t02, d, d0 (see sync file above)

dim(syncfile)
write.table(syncfile, paste("syncfile",args[4],args[5],args[6],args[8],"sync",sep="."), col.names = F, row.names = F, quote=F, sep="\t")

print("finished: creating sync file for Popoolationd")



