CCA <- function(file, Pfile, Sfile, groupFile, geneListSortFile = NULL, genome.build = c("Ch37","Ch38"), sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', plot = TRUE) 
{
  library(data.table) 
  library(RNMF)  
  mutdata = as.data.frame(fread(file))[, c(Hugo, chr, pos, ref, alt, sample.id)]
  ix = grepl('_', mutdata[,chr]) | grepl('M', mutdata[,chr])
  mut = mutdata[!ix,]
  mut[,chr] = sub('chr', '', mut[,chr])
  mut[,chr] = sub('Chr', '', mut[,chr])
  mut[,chr] = sub('CHR', '', mut[,chr])
  mut[,chr] = paste('chr', mut[,chr], sep='') 
  Nduplicate = mut[, c(Hugo, chr)]
  Nduplicate = unique(Nduplicate[order(Nduplicate[, Hugo]), ])
  Nduplicatenum = t(table(Nduplicate[, Hugo]))
  Nduplicate = colnames(Nduplicatenum)[which(Nduplicatenum>1)]
  index_ix = NULL
  for( i in 1:length(Nduplicate) )
  {
     ix = which(mut[, Hugo] == Nduplicate[i])
     if(length(ix)>0)
	 {
	    index_ix = c(index_ix, ix)
	 }
  }
  if(!is.null(index_ix))
  {
     mut = mut[-index_ix, ]
  }
  AnalCOSMICSigType = 'SBS'
  ID.mc.cores = 10
  ID.row83 = FALSE
  samMatrixsigsData = samMatrixSigs(Pfile, Sfile, AnalCOSMICSigType)
  geneMatrixNumData = geneMatrixNum(mut, geneListSortFile, 'MAF', AnalCOSMICSigType, genome.build, sample.id, chr, pos, ref, alt, Hugo, ID.mc.cores, ID.row83)
  geneList = geneMatrixNumData$geneList
  
  dir.create('samplesResults', showWarnings = FALSE, recursive = TRUE, mode = "0777")
  Samples = unlist(samMatrixsigsData$SampleID)
  n = length(Samples)
  if(is.null(groupFile))
  {  
    group = as.data.frame(cbind(Samples,'All'))
  }else{
    library(data.table)
    group = as.data.frame(fread(groupFile))
  }
  colnames(group) = c('SampleID', 'Group')
  ngs = unique(sort(group$Group))
  ngl = length(ngs)
  asnumeric = function(x)
  {
    x = as.numeric(as.character(x))
    return(x)
  }
  xy = function(x)
  {
    a = sum(x)
    if(a == 0)
    {
      a = 1
    }
    x = x/a
    return(x)
  }

  geneCumulativeContributionAbundance = NULL
  samMatrixDateRes = list()
  id_sam = 1
  geneListtmp = mut[, c(Hugo, chr)]
  data(knownGene)
  geneList1 = unique(merge(geneListtmp, knownGene, by=1:2, all.x = TRUE, sort = FALSE)) 
  colnames(geneList1) = c(Hugo, chr, 'Length')
  GinNs = intersect(geneList[,Hugo], geneList1[,Hugo])  
  for(i in 1:ngl)
  {
    dir.create(paste('samplesResults/image_',ngs[i],sep=""), showWarnings = FALSE, recursive = TRUE, mode = "0777")
    ix = group$Group == ngs[i]
    ngSamples = group$SampleID[ix]
    for(j in 1:length(ngSamples))
    {
	  cat('i:',i,'j:',j,'sample',ngSamples[j],'\n')
      samMatrixDateRes$SampleID[[id_sam]] = ngSamples[j]	
      dir.create(paste('samplesResults/', ngSamples[j], "/", AnalCOSMICSigType, sep=""), showWarnings = FALSE, recursive = TRUE, mode = "0777")
      index = which(samMatrixsigsData$SampleID == ngSamples[j])
      rho = samMatrixsigsData$P[[index]]
      rho = rho[, -1]
      rho = apply(rho, 2, asnumeric)
      index = which(geneMatrixNumData$SampleID == ngSamples[j])
	  if(length(index)>0)
	  {
		  G = geneMatrixNumData$samGeneMatrix[[index]]
		  G = G[, -1]
		  G = apply(G, 2, asnumeric)	  
          G = t(G)  # gene * type		  
	  }else{
		  G =  matrix(0, nrow(geneList1), nrow(samMatrixsigsData$P[[1]]))
		  rownames(G) = geneList1[, Hugo]
	  }
	  
      # 1
      GG = apply(G, 2, xy) 
      theta = GG%*%rho	
      theta = theta[GinNs, ]
	  tmp = t(apply(theta, 1, xy))
	  thetatable = cbind(t(t(rownames(tmp))), tmp)
	  colnames(thetatable)[1] = 'GeneName'
      write.table(thetatable, file=paste('samplesResults/', ngSamples[j], "/", AnalCOSMICSigType, "/", ngSamples[j], ".geneSigsResult.txt",sep=""), quote=F, col.names=T, row.names=F, sep="\t")
      samMatrixDateRes$theta[[id_sam]] = theta
	  
	  if(is.null(geneCumulativeContributionAbundance))
      {
        geneCumulativeContributionAbundance = theta
      }else{
        geneCumulativeContributionAbundance = geneCumulativeContributionAbundance + theta
      }

      if(plot)
      {
        theta = apply(theta, 1, xy)
        kn = nrow(theta)
        Nnameset = rownames(theta)
        for(sigi in 1:kn)
        {		
          xxda = theta[sigi, ]
		  imageMatrix = SMM(xxda)		 
          png(filename = paste('samplesResults/image_',ngs[i],"/",ngSamples[j],".",AnalCOSMICSigType,'.',Nnameset[sigi],'.png',sep=""),
              width = 448, height = 448)
          breaks.frequency <- seq(from=min(imageMatrix), to=max(imageMatrix), length.out=10)
          myColors <- colorRampPalette(c("grey98", "red"))
          image(1:nrow(imageMatrix), 1:ncol(imageMatrix), as.matrix(imageMatrix), breaks=breaks.frequency,
                col=myColors(length(breaks.frequency)-1), 
                axes = FALSE, cex = 1.5, xlab = "", ylab = "")
          dev.off()          
        }
      }
	  id_sam = id_sam + 1
    }
  }

  geneCumulativeContributionAbundancetable = cbind(t(t(rownames(geneCumulativeContributionAbundance))), geneCumulativeContributionAbundance)
  colnames(geneCumulativeContributionAbundancetable)[1] = 'GeneName'   
  write.table(geneCumulativeContributionAbundancetable, file=paste('samplesResults/', AnalCOSMICSigType, ".geneCumulativeContributionAbundance.txt",sep=""), quote=F, col.names=T, row.names=F, sep="\t")
}
