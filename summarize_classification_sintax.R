library(xlsx)

filenames <- list.files("reads-classified_sintax", pattern="*.txt", full.names=TRUE)
data <- lapply(filenames, read.table, header=F, sep="\t")

names <- gsub("reads-classified_sintax/", "", filenames)
names <- gsub(".sintax.txt", "", names)
names <- do.call('rbind', strsplit(as.character(names),'_',fixed=TRUE))[,1]
names(data) <- names

classification <- list()
for (sample in names(data)) {
  id <- data.frame(do.call('rbind', strsplit(as.character(data[[sample]]$V1),'_',fixed=TRUE)))[1]
  sequence <- data.frame(do.call('rbind', strsplit(as.character(data[[sample]]$V1),';',fixed=TRUE)))[1]
  size <- gsub(";", "", data[[sample]]$V1)
  size <- data.frame(do.call('rbind', strsplit(size,'=',fixed=TRUE)))[2]
  size$X2 <- as.numeric(size$X2)
  
  classif <- data.frame(do.call('rbind', strsplit(as.character(data[[sample]]$V2),',',fixed=TRUE)))
  k <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X1),'(',fixed=TRUE)))
  k$X1 <- gsub("k:", "", k$X1)
  k$X2 <- as.numeric(as.character(k$X2))
  p <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X2),'(',fixed=TRUE)))
  p$X1 <- gsub("p:", "", p$X1)
  p$X2 <- as.numeric(as.character(p$X2))
  c <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X3),'(',fixed=TRUE)))
  c$X1 <- gsub("c:", "", c$X1)
  c$X2 <- as.numeric(as.character(c$X2))
  o <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X4),'(',fixed=TRUE)))
  o$X1 <- gsub("o:", "", o$X1)
  o$X2 <- as.numeric(as.character(o$X2))
  f <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X5),'(',fixed=TRUE)))
  f$X1 <- gsub("f:", "", f$X1)
  f$X2 <- as.numeric(as.character(f$X2))
  g <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X6),'(',fixed=TRUE)))
  g$X1 <- gsub("g:", "", g$X1)
  g$X2 <- as.numeric(as.character(g$X2))
  s <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif$X7),'(',fixed=TRUE)))
  s$X1 <- gsub("s:", "", s$X1)
  s$X2 <- as.numeric(as.character(g$X2))
  
  classif <- cbind(id, sequence, k, p, c, o, f, g, s, size)
  colnames(classif)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)] <- c("sample", "sequence", "kingdom", "prob_kingdom", "division", "prob_division", "clade", "prob_clade", "order", "prob_order", "family", "prob_family", "genus", "prob_genus", "species", "prob_species", "size")
  
  rownames(classif) <- c()
  classification[[sample]] <- classif
}
rm(data) #clear memory

classification_all <- do.call("rbind", classification)
rm(classification) #clear memory
rownames(classification_all) <- c()

#filter out reads classified with probability <95% 
classification_all_filtered <- classification_all[classification_all$prob_species >= 0.95,] #& classification_all$size>=100

#aggregate Musa acuminata and Musa cultivar
classification_all_filtered[classification_all_filtered$species=="Musa hybrid cultivar_491892",]$species <- "Musa acuminata_4641"
#synonimize Eremobium lineare_1465652 and Eremobium aegyptiacum_368998
classification_all_filtered[classification_all_filtered$species=="Eremobium lineare_1465652",]$species <- "Eremobium aegyptiacum_368998"
#filter errorneous Alternanthera sp. XF30_1225511
classification_all_filtered <- classification_all_filtered[classification_all_filtered$species!="Alternanthera sp. XF30_1225511",]

#aggregate the reads after filtering
aggregated_all_filtered <- aggregate(size ~ sample+species, FUN = sum, classification_all_filtered)

#filter out reads with less than 100 copies per butterfly sample
aggregated_all_filtered <- aggregated_all_filtered[aggregated_all_filtered$size>100,]

#prepare data in a tabular form
aggregated_all_filtered$species <- as.factor(aggregated_all_filtered$species)
aggregated_all_filtered_table <- with(aggregated_all_filtered, tapply(size,list(species,sample),sum))
aggregated_all_filtered_table[is.na(aggregated_all_filtered_table)] <- 0
aggregated_all_filtered_table <- aggregated_all_filtered_table[,colSums(aggregated_all_filtered_table)!=0]

Total_samples <- rowSums(aggregated_all_filtered_table != 0)
Total_sequences <- rowSums(aggregated_all_filtered_table)

aggregated_all_filtered_table <- cbind(Total_sequences, Total_samples, aggregated_all_filtered_table)

#save as xlsx file
write.xlsx(aggregated_all_filtered_table, "aggregated_all_filtered_table.xlsx", row.names=T)
