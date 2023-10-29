library(SMimage)
Pfile = 'SBS.Pmatrix.txt'
Sfile = 'SBS.Smatrix.txt'
maffile = 'ESCC.exon.WES.maf-input-nonsilent.txt'
gfile = 'GroupfileNEW.txt'
genefile = 'genelistNew.txt'

CCA(file = maffile, Pfile = Pfile, Sfile = Sfile, groupFile = gfile, geneListSortFile = genefile, genome.build = "Ch37", sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', plot = TRUE) 
