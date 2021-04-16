##############################################################################
#Bacteria Network
#################
#####################
#Edit each step with your files!

Amazon_16S_esv_to_taxa_table <- read.csv("~/Google Drive/UCD/Research/2_Nfix_rate_15N2_DNA_SIP/16S/Amazon_16S_esv_to_taxa_table.csv", stringsAsFactors=FALSE, strip.white=TRUE)
Amazon_16S_esv_to_taxa_table$Class[which(Amazon_16S_esv_to_taxa_table$Class == "Nitrospira")] <- "Nitrospiria"
#View(Amazon_16S_esv_to_taxa_table)
t.esv <- t(rare_out)
write.csv(t.esv, file = "~/Google Drive/UCD/Research/2_Nfix_rate_15N2_DNA_SIP/16S/t.dada_seq.csv")

ra.taxa <- matrix(rep(-999), ncol = 31, nrow=3166)
ra.taxa <- as.data.frame(ra.taxa)
colnames(ra.taxa) <- c(colnames(t.esv)[1:24], colnames(Amazon_16S_esv_to_taxa_table)[2:8])
row.names(ra.taxa) <-   Amazon_16S_esv_to_taxa_table[,1]
for (p in seq(1:nrow(Amazon_16S_esv_to_taxa_table))) {
  matcher <- Amazon_16S_esv_to_taxa_table[p,1]
  index <- which(row.names(t.esv) == matcher)
  
  ra.taxa[p,1:24] <- t.esv[index,1:24]
  temp <- Amazon_16S_esv_to_taxa_table[p,2:8]
  temp[is.na(temp)] <- ""
  ra.taxa[p,25:31] <- temp
}
#View(ra.taxa)
ra.taxa$originalnumber <- seq(1:nrow(ra.taxa))


pasteName <- rep(0,nrow(ra.taxa))
for (p in seq(1:nrow(ra.taxa))) {
  if (ra.taxa$Species[p] =="") {
    if (ra.taxa$Genus[p] == "") {
      if (ra.taxa$Family[p] == "") {
        if (ra.taxa$Order[p] == "") {
          if (ra.taxa$Class[p] == "") {
            if (ra.taxa$Phylum[p] == "") {
              if (ra.taxa$Kingdom[p] == "") {
              } else {
                ra.taxa$Phylum[p] <- paste(ra.taxa$Kingdom[p], "sp.", sep = " ")
              }
            } else {
              ra.taxa$Class[p] <- paste(ra.taxa$Phylum[p], "sp.", sep = " ")
            }
          } else {
            ra.taxa$Order[p] <- paste(ra.taxa$Class[p], "sp.", sep = " ")
          }
        } else {
          ra.taxa$Family[p] <- paste(ra.taxa$Order[p], "sp.", sep = " ")
        }
      } else {
        ra.taxa$Genus[p] <- paste(ra.taxa$Family[p], "sp.", sep = " ")
      }
    } else {
      ra.taxa$Species[p] <- paste(ra.taxa$Genus[p], "sp.", sep = " ")
    }
  } 
  pasteName[p] <- paste(ra.taxa$Kingdom[p],ra.taxa$Phylum[p], ra.taxa$Class[p], ra.taxa$Order[p], ra.taxa$Family[p], ra.taxa$Genus[p], ra.taxa$Species[p], sep = " ")
  
}
ra.taxa$pasteName <- pasteName

for.aggregation <- data.frame(ra.taxa[,1:24])
row.names(for.aggregation) <- ra.taxa$originalnumber
Reduced.ra <- aggregate(. ~ ra.taxa$pasteName, for.aggregation, FUN="sum")
Reduced.t <- matrix(rep(-999), ncol = 7, nrow=nrow(Reduced.ra))
Reduced.t <- as.data.frame(Reduced.t)
colnames(Reduced.t) <- colnames(Amazon_16S_esv_to_taxa_table)[2:8]
taxa <- ra.taxa[,25:31]
row.names(taxa) <- seq(1:nrow(taxa))
for (g in seq(1:nrow(Reduced.ra))) {
  matcher <- Reduced.ra[g,1]
  index <- which(pasteName == matcher)
  temp <- taxa[index[1],]
  Reduced.t[g,] <- temp
}
Reduced <- cbind(Reduced.ra, Reduced.t)
Reduced <- Reduced[,-1]
nrow(Reduced)
ncol(Reduced)
#View(Reduced) #check taxa alignment
#Should be reduced to no repeated taxa
#New array is Reduced
write.table(Reduced, file = "~/Google Drive/UCD/Research/2_Nfix_rate_15N2_DNA_SIP/16S//Reduced_16_esv.csv")
backup_reduced_OTUs <- Reduced
#Now combine species with genus name and replace in existing matrix
New_spec_name <- -999.99
for(q in seq(1:nrow(Reduced))) {
  Spec <- Reduced[q,31]
  if (Spec != "" && grepl("sp.", Spec) == FALSE) {
    Gen <- Reduced[q,30]
    New_Spec <- paste(Gen,Spec, sep = " ")
    Spec <- New_Spec
  }
  New_spec_name <- c(New_spec_name, Spec)
}
New_spec_name <- New_spec_name[-1]
Reduced[,31] <- New_spec_name
#View(Reduced)
#Now discriminate based on distribution across at least 5 sites
Reduced_distrib <-  t(data.frame(rep(-99999,31)))
colnames(Reduced_distrib) <- colnames(Reduced)[1:31]
for (w in seq(1:nrow(Reduced))) {
  a <- Reduced[w,1:24]
  count = 0
  for (t in seq(1:length(a))) {
    if (a[t] > 0) {
      count = count + 1
    }
  }
  if (count > 2) {
    Reduced_distrib <- rbind(Reduced_distrib, Reduced[w,])
  }
  else {
    
  }
}
Reduced_distrib <- Reduced_distrib[-1,]
nrow(Reduced_distrib)
#Now Discriminate on overall abundance meeting a specified representation cutoff (0.01%)
########original run value
#Reduced_distrib_per_0.05 <- Reduced_distrib_per
###############
Totals <- rowSums(Reduced_distrib[,1:24])
All_total <- sum(Totals)
Percent <- 100*Totals/All_total
Reduced_distrib_per <-  t(data.frame(rep(-99999,31)))
colnames(Reduced_distrib_per) <- colnames(Reduced_distrib)[1:31]
for (w in seq(1:nrow(Reduced_distrib))) {
  if(Percent[w] > 0.009) {
    Reduced_distrib_per <- rbind(Reduced_distrib_per, Reduced_distrib[w,])
  }
  else {
    
  }
}
Reduced_distrib_per <- Reduced_distrib_per[-1,]
nrow(Reduced_distrib_per)
####MANUAL EDITS REQUIRED TO GET RID OF BAD NAMING. EXPORT AND IMPORT. THESE INCLUDED NITROSPIRA(G) AND TK10(O)
#write.table(Reduced_distrib_per, file = "~/Google Drive/OSU/Master's Project Data/Research-Masters etc/Statistics_Data_Analysis/Germany/Bacteria_Reduced_OTUs_fixed_taxonomy.csv")
#Reduced_distrib_per <- read.csv("~/Google Drive/OSU/Master's Project Data/Research-Masters etc/Statistics_Data_Analysis/Germany/Bacteria_Reduced_OTUs_fixed_taxonomy.csv", stringsAsFactors=FALSE)
#Then add a tail "sp" to groups where taxonomy ends, to distinguish from other groups with go further

Sp.add <- Reduced_distrib_per[,25:31]
for(g in seq(1:(ncol(Sp.add)-1))) {
  for(r in seq(1:nrow(Sp.add))) {
    # if (Sp.add[r,g+1] == "") {
    #  prev <- Sp.add[r,g]
    #  if (grepl(" sp.", prev)== FALSE && prev != "") {
    #    redo <- paste(prev, "sp.", sep= " ")
    #   Sp.add[r, g+1] <- redo
    #  }
    #}
    prev <- Sp.add[r,g]
    if (Sp.add[r,g+1] == prev && Sp.add[r,g+1] != "") {
      redo <- paste(".", Sp.add[r,g+1], ".", sep = "")
      Sp.add[r, g+1] <- redo
    }
  }
}
#View(Sp.add)
Reduced_distrib_per[,25:31] <- Sp.add 
write.table(Reduced_distrib_per, file = "~/Google Drive/UCD/Research/2_Nfix_rate_15N2_DNA_SIP/16S/Amazon_Reduced_further_edits_16S.csv")
###Check if there are names repeated across taxonomic levels. If there are, just add a '.' '.' to either side of the name

#########################################################
#Then quantify, test, and label
#Unique lists at each trophic level:
U_P <- Reduced_distrib_per[!duplicated(Reduced_distrib_per[,26]),26]
U_C <- Reduced_distrib_per[!duplicated(Reduced_distrib_per[,27]),27]
U_O <- Reduced_distrib_per[!duplicated(Reduced_distrib_per[,28]),28]
U_F <- Reduced_distrib_per[!duplicated(Reduced_distrib_per[,29]),29]
U_G <- Reduced_distrib_per[!duplicated(Reduced_distrib_per[,30]),30]
U_S <- Reduced_distrib_per[!duplicated(Reduced_distrib_per[,31]),31]
unique_tax <- list(P = U_P, C = U_C, O= U_O, Fa = U_F, G = U_G,S = U_S);  
attrib <- as.matrix(rbind(rep(999.9, 5),rep(999.9, 5)))
colnames(attrib) <- c("total", "pval","code", "label", "diff")
taxa <- c("taxa", "taxa")
for(g in seq(1:length(unique_tax))) {
  Un <- as.data.frame(unique_tax[g])
  for(q in seq(1:nrow(Un))) { #######EDIT nrow(Un)
    Un_comp <- as.character(Un[q,1])
    input <- (rep(0,24))
    col_ind <- g+25
    if(Un_comp != "") {
      for (b in seq(1:nrow(Reduced_distrib_per))) {
        if (Un_comp == Reduced_distrib_per[b,col_ind]) {
          #print("whooie") ######EDIT
          input <- input + Reduced_distrib_per[b,1:24]
        }
      }
      total <- sum(input)
      compare <- wilcox.test(as.numeric(input[1:12]), as.numeric(input[13:24]), paired = TRUE, exact = FALSE, conf.int = FALSE, warnings=FALSE)
      pval <- compare$p.value #adjust pvalues a the end
      diff <- mean(as.numeric(input[13:24])-as.numeric(input[1:12]))
      if (pval < 0.005 & diff < 0) { #code 5 means higher pre-harvest
        code <- 5
      } else if (pval >= 0.005 & pval < 0.05 & diff < 0) {
        code <- 4
      } else if (pval < 0.005 & diff > 0) { #code 2 means higher post-harvest
        code <- 2
      } else if (pval >= 0.005 & pval < 0.05 & diff > 0) {
        code <- 3
      } else {
        code <- 1
      }
      if (total > 100) {
        label <- 1
      } else if (code == 5) {
        label <- 1
      } else if (code == 2) {
        label <-1
      } else {
        label <- 0
      }
      line <- cbind(total, pval, code, label, diff)
      attrib <- rbind(attrib, line)
      taxa <- c(taxa, Un_comp)
      
    } else {
      
    }
  }
}
attrib <- attrib[-1:-2,]
taxa <- taxa[-1:-2]
adj_code <- -999.9
pval_adj <- p.adjust(attrib[,2], method = "BH")
diffy <- attrib[,5]
for (b in seq(1:length(pval_adj))) {
  pval <- pval_adj[b]
  if (pval < 0.005 & diffy[b] < 0) { #code 5 means forest
    code <- 5
  } 
  else if (pval >= 0.005 & pval < 0.05 & diffy[b] < 0) {
    code <- 4
  } 
  else if (pval < 0.005 & diffy[b] > 0) { #code 2 means higher pasture
    code <- 2
  } 
  else if (pval >= 0.005 & pval < 0.05 & diffy[b] > 0) {
    code <- 3
  } 
  else {
    code <- 1
  }
  
  
  adj_code <- cbind(adj_code, code)
  code <- c()
}
adj_code <- adj_code[-1]

adj_attrib <- cbind(attrib[,1], pval_adj, adj_code, attrib[,4], attrib[,5])
#####
##Add Bacteria to the data frames

Kingdom_sum <- aggregate(Reduced_distrib_per[,1:24],list(Reduced_distrib_per[,25]), FUN = "sum" )
View(Kingdom_sum)
Totals <- rowSums(Kingdom_sum[,2:25])
Arch_add <- c(Totals[1], 1, 1, 1,0)
Bact_add <- c(Totals[2], 1, 1, 1,0)
adj_attrib <- rbind(adj_attrib, Bact_add, Arch_add)

taxa <- c(taxa, "Bacteria", "Archaea")
##Code for naming
code_name <- taxa
for(g in seq(1:nrow(adj_attrib))) {
  if(adj_attrib[g,4] == 0) {
    code_name[g] <- ""
  }
}
attrib_out <- as.data.frame(adj_attrib[,1:5],taxa)
attrib_out <- attrib_out[,-4]
attrib_out[,5] <- row.names(attrib_out)
row.names(attrib_out) <- c()
attrib_out[,6] <- code_name
colnames(attrib_out) <- c("Totals", "pval","code","diff","Shared name","label")

nrow(attrib_out)
write.table(attrib_out, file = "~/Google Drive/UCD/Research/2_Nfix_rate_15N2_DNA_SIP/16S/Amazon_network_attrib_out_FINAL.csv")