\name{CCA}
\alias{CCA}
\alias{cca}
\alias{CA}
\alias{ca}

\title{
	Cumulative contribution abundance of genes.
}
\description{
	We present a tool to calculate the cumulative contribution abundance of genes in cancer research.
}
\usage{
	cumulativeCA(file = NULL, Pfile = NULL, Sfile = NULL, groupFile = NULL, geneListSortFile = NULL, genome.build = c("Ch37","Ch38"), sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', plot = TRUE) 	
}
\arguments{
  \item{file}{       
     A file with non-silent mutations obtained by annotation of mutation locations by Oncotator or other annotation software.
  }
  \item{Pfile}{       
     A mutational signature matrix file P with a format like the output result of RNMF software.
  }
   \item{Sfile}{       
     A abundance fractions matrix S with a format like the output result of RNMF software.
  }
  \item{groupFile}{       
     Group files with two columns: 'SampleID', 'Group'.
  }    
  \item{geneListSortFile}{ 
	Sorted gene list with two columns: 'Hugo_Symbol', 'Chromosome', such as "KDSR      chr18".
  }     
  \item{genome.build}{       
     The version of the reference genome used in the mutation set.
  }  
  \item{sample.id}{       
     Column name in the mutation file corresponding to the Sample ID.	
  }
  \item{chr}{       
     Column name in the mutation file corresponding to the chromosome.
  }  
  \item{pos}{       
     Column name in the mutation file corresponding to the mutation position.
  }  
  \item{ref}{       
     Column name in the mutation file corresponding to the reference base.
  }
  \item{alt}{       
     Column name in the mutation file corresponding to the alternate base.
  }
  \item{Hugo}{       
     Column name in the mutation file corresponding to the gene.
  }
  \item{plot}{       
    Whether to plot the image. 
  }  
}
\author{
	Zhenzhang Li(design), Wen Luo(implementation, testing)
}
\examples{
	CCA(file = maffile, Pfile = Pfile, Sfile = Sfile, groupFile = gfile, geneListSortFile = genefile, genome.build = "Ch37", sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', plot = TRUE)
}
