
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
    #print(Samp)
    SeqNo <- paste(split2[6], (unlist(strsplit(split2[7], ";")))[1], sep = ":")
    count <- as.numeric((unlist(strsplit(split2[7], "=")))[2])
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
    colnames(out.df) <- c("SampleName","seqID", "Count", "BlastID", "PercentID", "Cluster", "BlastTax", "Domain","King", "Order", "Class", "Family", "Genus", "Species")
  }
  
  return(out.df)
  close(con)
  close(con2)
  
} #called in IDunclassifiedSeqs

IDunclassifiedSeqs = function(ID_file, Blastfile, taxfile) {
  conA = file(ID_file, "r") #read in connection
  nomatch.df <- c()
  blast.frame <- BlOut_process(Blastfile, taxfile)
  while ( TRUE ) {
    line = readLines(conA, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    split_line <- unlist(strsplit(line, ":")) 
    SeqNo <- split_line[6:7]
    SeqNo2 <- unlist(strsplit(SeqNo[2], ";"))[1]
    SeqNoo <- paste(SeqNo[1], SeqNo2, sep=":")
    if (!(SeqNoo %in% blast.frame[,2])) { #this means this sequence entry had no blast hits!
      Samp_name <- gsub(">","",split_line[1])
      Count <- unlist(strsplit(SeqNo[2], "="))[2]
      Entree <- c(Samp_name, SeqNoo, Count, rep("No Match", 4), rep("Unclassified", 7))
      nomatch.df <- rbind(nomatch.df, Entree)
    }

  }
  if (nrow(nomatch.df >= 1)) {
    colnames(nomatch.df) <- c("SampleName", "seqID", "Count", "BlastID", "PercentID", "Cluster", "BlastTax", "Domain","King", "Order", "Class", "Family", "Genus", "Species")
    out.df <- rbind(blast.frame, nomatch.df)
  }
  return(out.df)
  close(conA)

}

#ID_file <- "~/Desktop/seq.nochim.fa"
#Blastfile <- "~/Desktop/out.blast.tsv"
#taxfile <- "~/Desktop/taxonomy.tsv"
#blast.frame2 <- BlOut_process(Blastfile,taxfile)

total.frame <- IDunclassifiedSeqs(ID_file,Blastfile,taxfile)
