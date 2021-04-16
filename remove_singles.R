remove_singles = function(otu_table, sfilt) {
  conA = file(otu_table, "r") #read in connection
  conB = file(sfilt, "a")
  while ( TRUE ) {
    line = readLines(conA, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(grepl('#',line)==TRUE) {
      writeLines(line, conB)
    } else {
      otu <- unlist(strsplit(line, '\t'))
      nummy <- sapply(otu[2:length(otu)], as.numeric)
      compr <- sum(nummy)
      if(compr > 1) {
        writeLines(line, conB)
      }
    }

    
  }
  close(conA)
  close(conB)
}

otu_table <- "~/all.otuout.fa"
sfilt <- "~/otu_table_no_single.tsv"
remove_singles(otu_table, sfilt)
 
