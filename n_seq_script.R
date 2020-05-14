#!/path/to/Rscript


#library to use getops
if (!requireNamespace("getopt", quietly = TRUE)){
  install.packages("getopt", repos='https://cloud.r-project.org')
}

library('getopt')

#help option as a function
help_opt = function() {
  cat("-v version\n-h help\n-i input file\n-n number of sequences we want
-t produce tree\n-d run diversity method\n-a aligned sequences file\n-p project\n")
}

#Version
ver = "2.0"

#getopt matrix
spec = matrix(c(
  'version', 'v', 0, 'double',
  'help', 'h', 0, 'character', 
  'input', 'i', 1, 'character',
  'num', 'n', 2, 'integer', 
  'tree', 't', 0, 'character',
  'cluster', 'c', 0, 'integer',
  'aligned', 'a', 0, 'character',
  'project', 'p', 1, 'character',
  'diversity', 'd', 0, 'integer'
), byrow=TRUE, ncol=4)
opt = getopt(spec)


#help option
if (!is.null(opt$help) ) {
  help_opt()
  q(status=1)
}


#version option
if (!is.null(opt$version) ) {
  cat("Version: ")
  cat(ver)
  cat("\n")
  q(status=1)
}

#libraries needed for sequence alignment
if (!requireNamespace("seqinr", quietly = TRUE)){
  install.packages("seqinr", repos='https://cloud.r-project.org')
}
if (!requireNamespace("DECIPHER", quietly = TRUE)){
  #library to install BioConductor
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='https://cloud.r-project.org')
  BiocManager::install()
  BiocManager::install("DECIPHER")
}
if (!requireNamespace("ape", quietly = TRUE)){
  install.packages("ape", repos='https://cloud.r-project.org')
}
if (!requireNamespace("factoextra", quietly = TRUE)){
  install.packages("factoextra", repos='https://cloud.r-project.org')
}

library("seqinr")
library("DECIPHER")
library("ape")
library("factoextra")


#input option
if (!is.null(opt$input) ){
  seqs = readDNAStringSet(opt$input)
  seqs = OrientNucleotides(seqs) #orient nucleotides in right direction if not already like that
  aligned_seqs = AlignSeqs(seqs) #perform alignment
}else{
  help_opt()
  q(status=1)
}

#num option
if (!is.null(opt$num) ){
  n = opt$num
}else{
  n = 1
}

#calculate distance matrix
dist = DistanceMatrix(aligned_seqs, type = "matrix", includeTerminalGaps = TRUE, 
                      penalizeGapGapMatches = FALSE, penalizeGapLetterMatches = FALSE)

#tree option
if (!is.null(opt$tree) ) {
  #Create tree by clustering sequences
  dend = IdClusters(dist, showPlot = TRUE, type = "dendrogram")
  #either convert to tree with write dendogram or save afa (aligned fasta)
  # saving the algined fasta is option two (save aligned file and use cod eon website to create tree)
  
  #Apply function to each node
  dend = dendrapply(dend, FUN = function(n) {
    if(is.leaf(n)) 
      attr(n, "label") <- 
        as.expression(substitute(italic(leaf),
                                 list(leaf=attr(n, "label"))))
    n
  })
  #write dendrogram to file then read in as tree
  if(!is.null(opt$project)){
    WriteDendrogram(dend, file = paste(opt$project, "_dend"))
    tree = read.tree(file = paste(opt$project, "_dend"))
  }else{
    WriteDendrogram(dend, file = "out_dend")
    tree = read.tree(file="out_dend") 
  }
 
  #Display tree
  if(!is.null(opt$project)){
    jpeg(filename = paste(opt$project, "_tree.jpg"), width = 500, height = 350)
  }else{
    jpeg(filename = "out_tree.jpg", width = 500, height = 350)
  }
  p = par(mar=c(1, 1, 1, 10),
           xpd=TRUE)
  plot(tree,
       yaxt="n")
  #last_plot.phylo = dend
  nodelabels(node = node.height(tree), col = "red", cex = 1.3)
  arrows(-0.1, 6, -0.2, 6,
         angle=90,
         length=0.05,
         code=3)
  text(-0.15, 6,
       "0.1",
       adj=c(0.5, -0.5))
  par(p)
  dev.off()
}

#diversity and number of sequences options
if (!is.null(opt$diversity) ) {
    final_set <- dist
    if(n == 1){
      cdists = rowSums(final_set)
      farthest <- which(cdists == max(cdists))[1]
      if(!is.null(opt$project)){
        write.table(farthest, file = paste(opt$project, "_diverse_sequences.txt"),  col.names = FALSE )
      }else{
        write.table(farthest, file = "out_diverse_sequences.txt", col.names = FALSE)
      }
    }else{
      while (nrow(final_set) > n) { 
        cdists = rowSums(final_set)
        closest <- which(cdists == min(cdists))[1]
        final_set <- final_set[-closest,-closest]
      }
      if(!is.null(opt$project)){
        write.table(row.names(final_set), file = paste(opt$project, "_diverse_sequences.txt"), sep = '\n', row.names = FALSE, col.names = FALSE)
      }else{
        write.table(row.names(final_set), file = "out_diverse_sequences.txt", sep = '\n', row.names = FALSE, col.names = FALSE)
      }
    }
}

if (!is.null(opt$cluster) ){
  clusters = hcut(dist, k = n)
  seq_list = c()
  for(i in 1:n){
    final_set <- dist[,colnames(dist) %in% names(clusters$cluster[clusters$cluster[] == i])]
    if(is.null(ncol(final_set))){
      seq_list <- append(seq_list, names(clusters$cluster[clusters$cluster[] == i]))
    }else{
      cdists <- colSums(final_set)
      farthest <- which(cdists == max(cdists))
      seq_list <- append(seq_list, names(farthest))
    }
  }
  
  if(!is.null(opt$project)){
    write.table(seq_list, file = paste(opt$project, "_diverse_sequences.txt"), sep = '\n', row.names = FALSE, col.names = FALSE)
  }else{
    write.table(seq_list, file = "out_diverse_sequences.txt", sep = '\n', row.names = FALSE, col.names = FALSE)
  }
}



#aligned option
if (!is.null(opt$aligned) ){
  if(!is.null(opt$project)){
    writeXStringSet(aligned_seqs, file = paste(opt$project, "_aligned_seqs")) 
  }else{
    writeXStringSet(aligned_seqs, file = "out_aligned_seqs") 
  }
}


#if nothing selected
if(is.null(opt$num) & is.null(opt$tree) & is.null(opt$div) & 
   is.null(opt$aligned) & is.null(opt$project)){
  help_opt()
}

