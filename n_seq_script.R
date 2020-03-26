#!/path/to/Rscript

#library to use getops
if (!requireNamespace("getopt", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("getopt")
}

library('getopt')

#help option as a function
help_opt = function() {
  cat("-v version\n-h help\n-i input file\n-n number of sequences we want
-t produce tree\n-d run diversity method\n-o output file\n")
}

#Version
ver = ".01"

#getopt matrix
spec = matrix(c(
  'version', 'v', 0, 'double',
  'help', 'h', 0, 'character', 
  'input', 'i', 2, 'character',
  'num', 'n', 2, 'integer', 
  'tree', 't', 0, 'character',
  'div', 'd', 2, 'integer',
  'output', 'o', 1, 'character'
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
  install.packages("BiocManager")
  BiocManager::install("seqinr")
}
if (!requireNamespace("DECIPHER", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("DECIPHER")
}
if (!requireNamespace("ape", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("ape")
}

library("seqinr")
library("DECIPHER")
library("ape")


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
  
  WriteDendrogram(dend, file = "dendrogram")
  tree = read.tree(file="dendrogram") #option one, save dendrogram to a file then read in the file as a tree.
  
  #Display tree
  jpeg(filename = "test_tree.jpg", width = 500, height = 350)
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
if (!is.null(opt$div) ) {
    final_set <- dist
    while (nrow(final_set) > n) {
      cdists = rowSums(final_set)
      closest <- which(cdists == min(cdists))[1]
      final_set <- final_set[-closest,-closest]
    }
    write.table(row.names(final_set), file = "diverse_sequences.txt", sep = '\n', row.names = FALSE, col.names = FALSE)
}


#output option
if (!is.null(opt$output) ){
  writeXStringSet(aligned_seqs, file= opt$output) #write to new file
  q(status=1)
}else{
  help_opt()
  q(status=1)
}


#if nothing selected
help_opt()

