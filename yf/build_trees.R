#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-a", "--alpha"), type="double", default=0.005, help="K-mer score threshold (logarithmic); default 0.005"),
  make_option(c("-c", "--cytoscape"), type="logical", default=FALSE, help="Create files for Cytoscape"),
  make_option(c("-f", "--full"), type="logical", default=TRUE, help="Full-length BCR mode; if V and J segments available"),
  make_option(c("-i", "--input"), help="Path to input file"),
  make_option(c("-o", "--output"), default="", help="Path to output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))


library(stringr)
library(dplyr)
library(fitdistrplus)
library(doParallel)
library(foreach)
library(igraph)
library(muscle)
library(Biostrings)
library(tidyr)

find_kmers <- function(string, k=5){
  n <- nchar(string) - k + 1
  kmers <- str_sub(string, 1:n, 1:n + k - 1)
  return(kmers)
}

interleave <- function(v1, v2){
  z <- c()
  for (i in 1:length(v1)){
    z <- c(z, v1[i], v2[i])
  }
  return(z)
}

merge_gaps <- function(mat, gaps, i){
  #count length of indel and merge gaps in an alignment matrix
  indel.num <- nrow(gaps[[i]])
  v1 <- mat[,i]
  v2 <- mat[,-i]
  if (indel.num > 0){
    for (x in 1:indel.num){
      v2 <- c( v2[0:(gaps[[i]][x,'start']-1)], gaps[[i]][x,'end']-gaps[[i]][x,'start']+1, 
               v2[(gaps[[i]][x,'end']+1):(length(v2)+1)] )
      v1 <- c( v1[0:(gaps[[i]][x,'start']-1)], '-', 
               v1[(gaps[[i]][x,'end']+1):(length(v1)+1)] )
      v2 <- v2[!is.na(v2)]
      v1 <- v1[1:length(v2)]
      gaps[[i]] <- gaps[[i]] - (gaps[[i]][x,'end']-gaps[[i]][x,'start'])
    }
  }
  matrix(c(v1,v2), ncol=2)
}

align2table <- function(alignment){
  #convert MultipleAlignment object to table with merged gaps
  gaps <- sapply(alignment, function(x) str_locate_all(x, "-+")[1])
  mat <- sapply(alignment, function(x) str_split(x, '')[[1]])
  for (i in 1:2){
    mat <- merge_gaps(mat, gaps, i)
  }
  data.frame(s1 = mat[,2], s2 = mat[,1])
}

mutations_weight <- function(mutations){
  prod(merge(data.frame(mutation = mutations), frequencies)$freq)
}

get_mutations_and_parent <- function(s1, s2){
  DNAset = DNAStringSet(c(s1, s2))
  alignment <- unmasked(muscle(DNAset, quiet = TRUE))
  aln.table <- align2table(alignment) %>%
    mutate(s1 = as.character(s1), s2 = as.character(s2))
  
  #if there are no indels
  if (!('-' %in% c(aln.table$s1, aln.table$s2))){
    ind <- which(aln.table$s1 != aln.table$s2)
    n1.aa <- c()
    n2.aa <- c()
    for (i in ind){
      codon_start = i - i%/%3 + 1
      codon1 = str_sub(s1, codon_start, codon_start + 2)
      codon2 = str_sub(s2, codon_start, codon_start + 2)
      aa1 = as.character(translate(DNAString(codon1)))
      aa2 = as.character(translate(DNAString(codon2)))
      n1.aa <- c(n1.aa, aa1)
      n2.aa <- c(n2.aa, aa2)
    }
    m <- aln.table[ind,] %>%
      mutate(n1 = paste(s1,s2,sep='>'), n2 = paste(s2,s1,sep='>'))
    m$ind <- ind
    m$n1.aa <- n1.aa
    m$n2.aa <- n2.aa
    
  } else {
    m <- aln.table %>% filter(s1 != s2) %>%
      mutate(n1 = ifelse(s1 == '-', paste0('D', s2), 
                         ifelse(s2 == '-', paste0('I', s1), paste(s1, s2, sep='>'))),
             n2 = ifelse(s2 == '-', paste0('D', s1), 
                         ifelse(s1 == '-', paste0('I', s2), paste(s2, s1, sep='>'))))
    m$n1.aa <- NA
    m$n2.aa <- NA
    m$ind <- NA
  }
  
  score1 <- mutations_weight(m$n1)
  score2 <- mutations_weight(m$n2)
  parent = ifelse(score1 > score2, 'n1', 'n2')
  child = ifelse(parent=='n1', 'n2', 'n1')
  m$mut.nt <- m[[child]]
  list(parent = parent, mutations = dplyr::select(m, -s1, -s2, -n1, -n2))
}

pair_compare <- function(i, full, alpha){
  # find all clonotypes related to clonotype number i
  n1.lst <- c()
  n2.lst <- c()
  kmer.score.lst <- c()
  mutations.num.lst <- c()
  cdr3.muts.lst <- list()
  n1.aa.lst <- list()
  n2.aa.lst <- list()
  position.lst <- list()
  mut.type.lst <- list()
  
  df <- filter(.df, v == .df$v[i], j == .df$j[i])
  
  for (j in df$clone[df$clone > i]){
    #cat(i,' - ',j,'\n')
    shared.kmers = intersect(kmers[[i]], kmers[[j]])
    shared.kmer.inf = sum(unlist(lapply(shared.kmers, function(x) kmer.df$neg_ln[kmer.df$kmer == x])))
    inf1 = sum(unlist(lapply(kmers[[i]], function(x) kmer.df$neg_ln[kmer.df$kmer == x])))
    inf2 = sum(unlist(lapply(kmers[[j]], function(x) kmer.df$neg_ln[kmer.df$kmer == x])))
    m = shared.kmer.inf/(inf1*inf2)
    #p.value = 1 - pgamma(shared.kmer.inf/(inf1*inf2), shape = gamma.params$estimate['shape'], rate = gamma.params$estimate['rate'])
      
    if (log(m) > alpha){
      if (full){ #full-length mode
        mut.i = df[df$clone == i,]$all.mutations
        mut.j = df[df$clone == j,]$all.mutations
        shared.muts = intersect(mut.i, mut.j)
          
        if (length(shared.muts) == length(mut.i) | length(shared.muts) == length(mut.j)){
          n1 = ifelse(length(shared.muts) == length(mut.i),i, j)
          n2 = ifelse(n1 == i, j, i)
          cdr3.muts <- list(get_mutations(df[df$clone == n1,]$cdr3nt, df[df$clone == n2,]$cdr3nt)$n2)
          muts.between = list(setdiff(df[df$clone == n2,]$all.mutations, df[df$clone == n1,]$all.mutations))
            
          cdr3.muts.lst[[paste(n1,n2,sep='_')]] <- cdr3.muts
          non.cdr3.muts.lst[[paste(n1,n2,sep='_')]] <- muts.between
          n1.lst <- c(n1.lst, n1)
          n2.lst <- c(n2.lst, n2)
          kmer.score.lst <- c(kmer.score.lst, m)
          mutations.num.lst <- c(mutations.num.lst, mutations.num = 
                                   length(muts.between[[1]]) + length(cdr3.muts[[1]]))
        }
      }
        
      else {
        dp = get_mutations_and_parent( df[df$clone == i,]$cdr3nt, df[df$clone == j,]$cdr3nt )
          
        n1 = ifelse(dp$parent == 'n1', i, j)
        n2 = ifelse(dp$parent == 'n1', j, i)
        n1.lst <- c(n1.lst, n1)
        n2.lst <- c(n2.lst, n2)
        kmer.score.lst <- c(kmer.score.lst, m)
        mutations.num.lst <- c(mutations.num.lst, mutations.num = length(dp$mutations$mut.nt))
        name = paste(n1,n2,sep='_')
          
        if (length(dp$mutations$mut.nt) == 0){
          cdr3.muts.lst[[name]] <- NA
          n1.aa.lst[[name]] <- NA
          n2.aa.lst[[name]] <- NA
          position.lst[[name]] <- NA
          mut.type.lst[[name]] <- NA
        } else {
          cdr3.muts.lst[[name]] <- dp$mutations$mut.nt
          n1.aa.lst[[name]] <- dp$mutations$n1.aa
          n2.aa.lst[[name]] <- dp$mutations$n2.aa
          position.lst[[name]] <- dp$mutations$ind
          mut.type.lst[[name]] <- ifelse(dp$mutations$n1.aa == dp$mutations$n2.aa, 'S', 'R')
        }
      }
    }
  }
  
  if (full){
    data.frame(n1 = n1.lst, n2 = n2.lst, kmer.score = kmer.score.lst, 
               mut.num = mutations.num.lst, cdr3.muts = I(cdr3.muts.lst), 
               non.cdr3.muts = I(non.cdr3.muts.lst))
               #n1.aa = I(n1.aa.lst), 
               #n2.aa = I(n2.aa.lst), position = I(position.lst),
               #mut.type = I(mut.type.lst))
  }
  else{
    data.frame(n1 = n1.lst, n2 = n2.lst, kmer.score = kmer.score.lst, 
               mut.num = mutations.num.lst, cdr3.muts = I(cdr3.muts.lst), 
               n1.aa = I(n1.aa.lst), n2.aa = I(n2.aa.lst), position = I(position.lst),
               mut.type = I(mut.type.lst))
  }
}

get_arbor_edge <- function(node){
  old <- filter(pairs, n2 == node)
  old[which.min(old$mut.num),]
}

clone_info <- function(tree, full){
  
  clone <- data.frame(root = numeric(), cdr3aa = character(), freq = double(), freq.sum = double(), 
                      leaves = integer(), nodes = integer(), root.mut = integer(), diameter = integer(),
                      mean.mut = double(), total.mut = integer(), branching = double(),
                      mean.degree = double())
  
  if (length(V(tree)) > 1){
    edges <- get.edgelist(tree)
    root <- setdiff(edges[,1], edges[,2])
    root.i = as.integer(root)
    cdr3aa <- df[root.i,]$cdr3aa
    freq <- df[root.i,]$freq
    freq.sum = sum(df[as.numeric(names(V(tree))),]$freq)
    
    leaves <- setdiff(edges[,2], edges[,1])
    mut.from.root <- distances(tree, root)
    
    if(full){
      root.mut <- length(df$all.mutations[[root.i]])
    } else{root.mut <- 0}
    
    diameter = max(mut.from.root)
    mean.mut = mean(mut.from.root)
    total.mut = sum(E(tree)$weight)
    
    branching = length(leaves)/mean.mut
    mean.degree = mean(degree(tree))
    
    clone <- rbind(clone, list(root = root, cdr3aa = cdr3aa, freq = freq, freq.sum = freq.sum, 
                               leaves = length(leaves), nodes = length(V(tree)), root.mut = root.mut,
                               diameter = diameter, mean.mut = mean.mut, total.mut = total.mut, 
                               branching = branching, mean.degree = mean.degree))
  }
  return(clone)
}

make_cytoscape_files <- function(sample, edge.table){
  
  # write net file
  net <- data.frame(from = edge.table$n1, to = edge.table$n2, interaction = rep('shm', nrow(edge.table)))
  write.table(net, file=paste0('cytoscape/', sample, '.', tissue, '.net.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  # write edge file
  edge <- data.frame(edge = paste(edge.table$n1, '(shm)', edge.table$n2, sep=' '), 
                     mut.num = edge.table$mut.num,
                     cdr3.muts = edge.table$cdr3.muts)
  write.table(edge, file=paste0('cytoscape/', sample, '.', tissue, '.edge.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  # write node file
  node <- df %>%
    filter(sample == sample, clone %in% all.nodes) %>%
    dplyr::select(node = clone, cdr3aa, v, j, freq, isotype)
  write.table(node, file=paste0('cytoscape/', sample, '.', tissue, '.node.txt'), sep='\t', row.names=FALSE, quote=FALSE)
}


#import frequencies of k-mers and mutations
load('~/yf/mut_frequencies.rda')
load('~/yf/kmer_frequencies.rda')

kmer.df$neg_ln <- -log(kmer.df$freq)

#prepare dataset
df <- read.table(opt$i, header=T, sep="\t")
df$v <- sapply(df$v, function(x) str_split(x, ',')[[1]][1])
df$j <- sapply(df$j, function(x) str_split(x, ',')[[1]][1])
#df$v <- str_sub(str_extract(df$v, '(.+)\\*'), 1, -2)
#df$j <- str_sub(str_extract(df$j, '(.+)\\*'), 1, -2)

if(opt$f){
  df <- df %>% 
    mutate(all.mutations = paste(mutations.nt.FR1, mutations.nt.CDR1, mutations.nt.FR2, mutations.nt.CDR2,
                                          mutations.nt.FR3, mutations.nt.CDR3, mutations.nt.FR4, sep=','),
           cdr3nt = as.character(cdr3nt)) %>%
    mutate(all.mutations = str_extract_all(all.mutations, '([\\w\\d>:]+)'),
           ndn = str_sub(cdr3nt, pmax(0, v.end.in.cdr3-4), pmin(j.start.in.cdr3+4, nchar(cdr3nt))))
  }

cl <- makeCluster(detectCores())
registerDoParallel(cl)

shm <- data.frame()
clones <- data.frame()

for (s in unique(df$sample)){
  
  .df <- filter(df, sample == s) %>%
    mutate(cdr3nt = as.character(cdr3nt))
  
  # list all 5-mers
  if(opt$f){
    kmers <- foreach(i = .df$ndn, .packages='stringr') %dopar% find_kmers(i)
  } else{
    kmers <- foreach(i = .df$cdr3nt, .packages='stringr') %dopar% find_kmers(i)
  }
  
  # get and filter edges
  pairs <- foreach(x = .df$clone, .combine='rbind', .packages = c('dplyr', 'stringr', 'Biostrings', 'muscle')) %dopar% pair_compare(x, full = opt$f, alpha = opt$a)
  final.pairs <- foreach(x = .df$clone, .combine='rbind', .packages = c('dplyr')) %dopar% get_arbor_edge(x)
  
  # graph analysis
  g <- graph( edges=as.character(interleave(final.pairs$n1, final.pairs$n2)) )
  g <- set.edge.attribute(g, 'weight', value = final.pairs$mut.num)
  components <- decompose.graph(g)
  .clones <- foreach(x = components, .combine='rbind', .packages = c('igraph')) %dopar% clone_info(x, full = opt$f)
  .clones$single <- FALSE
  .clones <- mutate(.clones, root = as.integer(root), cdr3aa = as.character(cdr3aa))
  
  if (opt$f){
    singletons <- .df[-as.numeric(names(V(g))), ] %>% dplyr::select(root = clone, cdr3aa, freq, all.mutations) %>%
      mutate(freq.sum = freq, leaves = 0, nodes = 1, diameter = 0,
             mean.mut = 0, total.mut = 0, branching = 0, mean.degree = 0, single = TRUE)
    singletons$root.mut <- sapply(singletons$all.mutations, length)
    singletons <- dplyr::select(singletons, -all.mutations)
  } else{
    singletons <- .df[-as.numeric(names(V(g))), ] %>% dplyr::select(root = clone, cdr3aa, freq) %>%
      mutate(freq.sum = freq, leaves = 0, nodes = 1, diameter = 0,
             mean.mut = 0, total.mut = 0, branching = 0, mean.degree = 0, single = TRUE)
    singletons$root.mut <- 0
  }
  
  .clones <- rbind(.clones, singletons)
  .clones$sample <- s
  clones <- rbind(clones, .clones)
  
  #shm - not ready yet
  if (opt$f){
    .shm <- dplyr::select(final.pairs, n2, non.cdr3.muts, v, j)
    .shm$muts <- sapply(.new.shm$non.cdr3.muts, unlist)
  } else{
    .shm <- data.frame()
  }
  .shm <- .shm %>% unnest(muts = muts) %>%
    mutate(pos = str_extract(muts, '\\d+'), from=str_sub(str_extract(muts, '.>'),1,1),
           to = str_sub(muts,-1,-1), sample = s) %>%
    dplyr::select(-muts)
  
  .cdr3.shm <- final.pairs[sapply(final.pairs$cdr3.muts, function(x) length(x[[1]])) > 0,]
  .cdr3.shm <- mutate(.cdr3.shm, sample=s, clonotype.id=n2, pos.nt=NA, segment.name=NA,
         region='CDR3', mut.type = unlist(.cdr3.shm$mut.type), position = unlist(.cdr3.shm$position),
             cdr3.mut = unlist(.cdr3.shm$cdr3.muts), from.aa = unlist(.cdr3.shm$n1.aa),
             to.aa = unlist(.cdr3.shm$n2.aa)) %>%
    filter(!is.na(to.aa))
  
  .shm$sample <- s
  shm <- rbind(cdr3.shm, .cdr3.shm)
  
  if (opt$c){
    make_cytoscape_files()
  }
}

stopCluster(cl)

