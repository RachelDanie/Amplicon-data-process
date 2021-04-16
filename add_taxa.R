#!/usr/bin/env Rscript

args= commandArgs (trailingOnly=TRUE)

BlOut_process = function(Blastfile, taxfile) {
  con = file(Blastfile, "r")
  con2 = file(taxfile, "r")#read in connection
  out.df <- c()
  TaxVect = readLines(con2)
  Bl_ID = c()
  Tax_line = c()
  for (p in seq(1:length(TaxVect))) {
    line <- unlist(strsplit(TaxVect[p], "\t"))
    Bl_ID <- c(Bl_ID,line[1])
    Tax_line <- c(Tax_line, line[4] )
  }
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) { 
      break
    }
    split1 <- unlist(strsplit(line, "\t")) 
    pID <- as.numeric(split1[4]) #Percent Identity
    split2 <- unlist(strsplit(split1[1], ":"))
    Samp <- split2[1]
    SeqNo <- paste(split2[6], (unlist(strsplit(split2[7], ";")))[1], sep = ":")
    count <- as.numeric((unlist(strsplit(split2[7], "size=")))[2])
    split3 <- unlist(strsplit(split1[2], ";"))
    BlastID <- split3[1]
    Cluster <- gsub("_", " ", split3[2])
    Spec_check <- split3[3]
    if (pID < 75) {
      taxLevel <- 4
    } else if (pID >= 75 && pID <= 88.1) {
      taxLevel = 5
    } else if (pID > 88.1 && pID <= 91.9) {
      taxLevel = 6
    } else {
      taxLevel = 7
    }
    
    if (BlastID %in% Bl_ID) {
      taxline <- unlist(strsplit(Tax_line[match(BlastID, Bl_ID)], ";"))
      Taxout <- taxline[1:taxLevel]
      if (length(Taxout) < 7) {
        Taxout <- c(Taxout, rep("Unclassified", 7-length(Taxout)))
      }
      
    }
    
    Entry <- c(Samp, SeqNo, count, BlastID, pID, Cluster, Spec_check, Taxout)
    out.df <- rbind(out.df, Entry)
    colnames(out.df) <- c("SampleName","seqID", "Count", "BlastID", "PercentID", "Cluster", "BlastTax", "King", "Phylum", "Order", "Class", "Family", "Genus", "Species")
  }
  
  return(out.df)
  close(con)
  close(con2)
  
} 


add_taxa = function(otu_table_filt, Blastfile, taxfile, out_final) {
  conA = file(otu_table_filt, "r") #read in connection
  conB = file(out_final, "a")
  blast.frame <- BlOut_process(Blastfile, taxfile)
  while ( TRUE ) {
    line = readLines(conA, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(grepl('#',line)==TRUE) {
      newLine <- paste(line, "seqID", "Count", "PercentID", "BlastTax", "Cluster", "King", "Phylumn", "Order", "Class", "Family", "Genus", "Species", sep="\t")
      writeLines(newLine, conB)
    } else {
      otu <- unlist(strsplit(line, '\t'))
      md <- unlist(strsplit(otu[1],':'))
      SequenceID <- paste(md[6], md[7], sep=":")
      Sampname <- md[1]
      if((SequenceID %in% blast.frame[,2])) { #If there was a BLAST HIT
        entr <- which(blast.frame[,2] == SequenceID)
        newLine <- paste(line, blast.frame[entr,2], blast.frame[entr,3], blast.frame[entr,5], blast.frame[entr,7], blast.frame[entr,6], blast.frame[entr,8], blast.frame[entr,9], blast.frame[entr,10], blast.frame[entr,11], blast.frame[entr,12], blast.frame[entr,13],blast.frame[entr,14], sep="\t")
        writeLines(newLine, conB)
      } else {
	cnt <- sum(as.numeric(unlist(strsplit(line, "\t"))[2:length(unlist(strsplit(line, "\t")))])) 
        newLine <- paste(line, "NA", cnt, "NA", "NA","Unclassified", "Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified",sep="\t")
        writeLines(newLine, conB)
      }
    }
    
    
  }
  close(conA)
  close(conB)
}

add_taxa(args[1], args[2], args[3], args[4])
