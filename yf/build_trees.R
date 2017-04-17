# to do: add positions for cdr3 mutations

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

mutations_weight <- function(mutations){
  prod(merge(data.frame(mutation = mutations), frequencies)$freq)
}

merge_gaps <- function(mat, gaps, i){
  #count length of indel and merge gaps in an alignment matrix
  indel.num <- nrow(gaps[[i]])
  v1 <- mat[,i]
  v2 <- mat[,-i]
  if (indel.num > 0){
    for (x in 1:indel.num){
      v2 <- c( v2[0:(gaps[[i]][x,'start']-1)], gaps[[i]][x,'end']-gaps[[i]][x,'start']+1, 
               v2[(gaps[[i]][x,'end']+1):length(v2)] )
      v1 <- c( v1[0:(gaps[[i]][x,'start']-1)], '-', v1[(gaps[[i]][x,'end']+1):length(v1)] )
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

get_mutations <- function(s1, s2){
  DNAset = DNAStringSet(c(s1, s2))
  alignment <- unmasked(muscle(DNAset, quiet = TRUE))
  aln.table <- align2table(alignment)
  
  aln.table %>% mutate(s1 = as.character(s1), s2 = as.character(s2)) %>%
    filter(s1 != s2) %>%
    mutate(n1 = ifelse(s1 == '-', paste0('D', s2), 
                       ifelse(s2 == '-', paste0('I', s1), paste(s1, s2, sep='>'))),
           n2 = ifelse(s2 == '-', paste0('D', s1), 
                       ifelse(s1 == '-', paste0('I', s2), paste(s2, s1, sep='>')))) %>%
    dplyr::select(-s1, -s2)
}

define_parent <- function(mut.table){
  score1 <- mutations_weight(mut.table$n1)
  score2 <- mutations_weight(mut.table$n2)
  
  parent = ifelse(score1 > score2, 'n1', 'n2')
  list(parent = parent, mutations=dplyr::select(cdr3.muts, -get(parent))[,1])
}

pair_compare <- function(i, gm = 'yes', alpha = 1e-08){
  # find all clonotypes related to clonotype number i
  cdr3.muts.lst <- list()
  non.cdr3.muts.lst <- list()
  n1.lst <- c()
  n2.lst <- c()
  kmer.p.value.lst <- c()
  mutations.num.lst <- c()
  
  for (j in (i+1):nrow(df)){
    #cat(i,' - ',j,'\n')
    if (df$v[i] == df$v[j] & df$j[i] == df$j[j]){
      shared.kmers = intersect(kmers[[i]], kmers[[j]])
      shared.kmer.inf = sum(unlist(lapply(shared.kmers, function(x) kmer.df$neg_ln[kmer.df$kmer == x])))
      inf1 = sum(unlist(lapply(kmers[[i]], function(x) kmer.df$neg_ln[kmer.df$kmer == x])))
      inf2 = sum(unlist(lapply(kmers[[j]], function(x) kmer.df$neg_ln[kmer.df$kmer == x])))
      p.value = 1 - pgamma(shared.kmer.inf/(inf1*inf2), shape = gamma.params$estimate['shape'], rate = gamma.params$estimate['rate'])
      
      if (p.value < alpha){
        
        if (gm == 'yes'){ #if germline mutations are available
          mut.i = df$all.mutations[[i]]
          mut.j = df$all.mutations[[j]]
          shared.muts = intersect(mut.i, mut.j)
          
          if (length(shared.muts) == length(mut.i) | length(shared.muts) == length(mut.j)){
            n1 = ifelse(length(shared.muts) == length(mut.i),i, j)
            n2 = ifelse(n1 == i, j, i)
            cdr3.muts <- list(get_mutations(df$cdr3nt[n1], df$cdr3nt[n2])$n2)
            muts.between = list(setdiff(df$all.mutations[[n2]], df$all.mutations[[n1]]))
            
            cdr3.muts.lst[[paste(n1,n2,sep='_')]] <- cdr3.muts
            non.cdr3.muts.lst[[paste(n1,n2,sep='_')]] <- muts.between
            n1.lst <- c(n1.lst, n1)
            n2.lst <- c(n2.lst, n2)
            kmer.p.value.lst <- c(kmer.p.value.lst, p.value)
            mutations.num.lst <- c(mutations.num.lst, mutations.num = 
                                     length(muts.between[[1]]) + length(cdr3.muts[[1]]))
          }
        }
        
        else {
          dp = define_parent(cdr3.muts)
          cdr3.muts <- dplyr::select(cdr3.muts, -get(dp$parent))[,1]
          row <- list(n1 = ifelse(dp$parent == 'n1', i, j), n2 = ifelse(dp$parent == 'n1', j, i),
                      kmer.p.value = p.value, mutations.num = dp$mutations.num)
          
          cdr3.muts.lst[[paste(n1,n2,sep='_')]] <- cdr3.muts
          n1.lst <- c(n1.lst, n1)
          n2.lst <- c(n2.lst, n2)
          kmer.p.value.lst <- c(kmer.p.value.lst, p.value)
          mutations.num.lst <- c(mutations.num.lst, mutations.num = length(dp$mutations))
          
        }
      }
    }
  }
  if (gm == 'yes'){
    data.frame(cdr3.muts = I(cdr3.muts.lst), non.cdr3.muts = I(non.cdr3.muts.lst),n1 = n1.lst,
               n2 = n2.lst, kmer.p.value = kmer.p.value.lst, mut.num = mutations.num.lst)
  }
  else{
    data.frame(cdr3.muts = I(cdr3.muts.lst), n1 = n1.lst,
               n2 = n2.lst, kmer.p.value = kmer.p.value.lst, mut.num = mutations.num.lst)
  }
}

get_arbor_edge <- function(node){
  old <- filter(pairs, n2 == node)
  old[which.min(old$mut.num),]
}

clone_info <- function(tree){
  clone <- data.frame(root = numeric(), cdr3nt = character(), freq = double(), freq.sum = double(), 
                      leaves = integer(), nodes = integer(), root.mut = integer(), diameter = integer(),
                      mean.mut = double(), total.mut = integer(), branching = double(),
                      mean.degree = double())
  
  if (length(V(tree)) > 1){
    edges <- get.edgelist(tree)
    root <- setdiff(edges[,1], edges[,2])
    root.i = as.integer(root)
    cdr3nt <- df[root.i,]$cdr3nt
    freq <- df[root.i,]$freq
    freq.sum = sum(df[as.numeric(names(V(tree))),]$freq)
    
    leaves <- setdiff(edges[,2], edges[,1])
    mut.from.root <- distances(tree, root)
    root.mut <- length(df$all.mutations[[root.i]])
    
    diameter = max(mut.from.root)
    mean.mut = mean(mut.from.root)
    total.mut = sum(E(tree)$weight)
    
    branching = length(leaves)/mean.mut
    mean.degree = mean(degree(tree))
    
    clone <- rbind(clone, list(root = root, cdr3nt = cdr3nt, freq = freq, freq.sum = freq.sum, 
                               leaves = length(leaves), nodes = length(V(tree)), root.mut = root.mut, 
                               diameter = diameter, mean.mut = mean.mut, total.mut = total.mut, 
                               branching = branching, mean.degree = mean.degree))
  }
  return(clone)
}

make_cytoscape_files <- function(edge.table){
  #write net file
  net <- data.frame(from = edge.table$n1, to = edge.table$n2, interaction = rep('shm', nrow(edge.table)))
  write.table(net, file=paste0('trees/', sample, '.net.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  # write edge file
  edge <- data.frame(edge = paste(edge.table$n1, '(shm)', edge.table$n2, sep=' '), mutations = edge.table$mutations)
  write.table(edge, file=paste0('trees/', sample, '.edge.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  # write node file
  all.nodes <- unique(c(edge.table$n1, edge.table$n2))
  node <- df[all.nodes,] %>% dplyr::select(cdr3aa, v, j, freq)
  node$node <- all.nodes
  write.table(node, file=paste0('trees/', sample, '.node.txt'), sep='\t', row.names=FALSE, quote=FALSE)
}




#import frequencies of k-mers and mutations and fit_gamma
load('~/yf/mut_frequencies.rda')
load('~/yf/kmer_frequencies.rda')
load('~/yf/gamma_params.rda')

kmer.df$neg_ln <- -log(kmer.df$freq)

old <- paste0('yf_old_RNA/', c("Abdulain", "Ilgen", "Mamaev", "Smirnov", "Vlasov"))
young <- paste0('yf_young_RNA/', c("Antipyat", "Epifancev", "Hadjibekov", "Koshkin", "Kovalchuk"))

cl <- makeCluster(detectCores())
registerDoParallel(cl)

new.shm <- data.frame()
cdr3.shm <- data.frame()

for (sample in c(old, young)){
  #prepare dataset
  df <- read.table(paste0('data/', sample, '.clones.txt'), header=T, sep="\t")
  df <- df %>% mutate(all.mutations=paste(mutations.nt.FR1, mutations.nt.CDR1, mutations.nt.FR2, mutations.nt.CDR2,
                                          mutations.nt.FR3, mutations.nt.CDR3, mutations.nt.FR4, sep=','),
                      cdr3nt = as.character(cdr3nt)) %>%
    mutate(all.mutations = str_extract_all(all.mutations, '([\\w\\d>:]+)'),
           ndn = str_sub(cdr3nt, pmax(0, v.end.in.cdr3-4), pmin(j.start.in.cdr3+4, nchar(cdr3nt))))
  #df <- df %>% group_by(cdr3nt, v.end.in.cdr3, j.start.in.cdr3, v, j, all.mutations)  %>% summarise(sample = 'raji')
  df$v <- str_sub(str_extract(df$v, '(.+)\\*'), 1, -2)
  df$j <- str_sub(str_extract(df$j, '(.+)\\*'), 1, -2)
  
  # list all 5-mers
  kmers <- foreach(i = df$ndn, .packages='stringr') %dopar% find_kmers(i)
  
  # get and filter edges
  pairs <- foreach(x = 1:(nrow(df)-1), .combine='rbind', .packages = c('dplyr', 'stringr', 'Biostrings', 'muscle')) %dopar% pair_compare(x)
  final.pairs <- foreach(x = 1:nrow(df), .combine='rbind', .packages = c('dplyr')) %dopar% get_arbor_edge(x)
  
  #graph analysis
  g <- graph( edges=as.character(interleave(final.pairs$n1, final.pairs$n2)) )
  g <- set.edge.attribute(g, 'weight', value = final.pairs$mut.num)
  components <- decompose.graph(g)
  clones <- foreach(x = components, .combine='rbind', .packages = c('igraph')) %dopar% clone_info(x)
  clones$single <- FALSE
  clones <- mutate(clones, root = as.integer(root), cdr3nt = as.character(cdr3nt))
  
  singletons <- df[-as.numeric(names(V(g))), ] %>% dplyr::select(cdr3nt, freq, all.mutations) %>%
    mutate(freq.sum = freq, leaves = 0, nodes = 1, diameter = 0,
           mean.mut = 0, total.mut = 0, branching = 0, mean.degree = 0, single = TRUE)
  singletons$root.mut <- sapply(singletons$all.mutations, length)
  singletons$root <- setdiff(1:nrow(df), as.numeric(names(V(g))))
  singletons <- dplyr::select(singletons, -all.mutations)
  
  clones <- rbind(clones, singletons)
  
  write.table(clones, file=paste0('trees/stat/', sample, '.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  #filter shm table
  type = str_extract(sample, '.NA')
  patient = str_sub(str_extract(sample, '\\/(.+)'), 2)
  proj = str_sub(str_extract(sample, 'yf_[^_]+_'), 4, -2)
    
  .new.shm <- dplyr::select(final.pairs, n2, non.cdr3.muts)
  .new.shm$pos <- sapply(.new.shm$non.cdr3.muts, function(x) str_extract(unlist(x), '\\d+'))
  .new.shm <- .new.shm %>% unnest(pos = pos) %>% mutate(type = type, sample = patient)
  new.shm <- rbind(new.shm, .new.shm)
  
  .cdr3.shm <- final.pairs[sapply(final.pairs$cdr3.muts, function(x) length(x[[1]])) > 0,]
  .cdr3.shm <- mutate(.cdr3.shm, proj=proj, sample=sample, type=type, replica=1, allele.rate=0,
                      allele.clonotype.share=0, clonotype.id=n2, pos.nt=NA, segment.name=NA,
                      region='CDR3')
  cdr3.shm <- rbind(cdr3.shm, .cdr3.shm)
}

stopCluster(cl)

#preprocessed shm table
load("sp.Rda")
shm <- df

new.shm <- transmute(new.shm, clonotype.id=n2-1, pos.nt=as.integer(pos), sample=sample, type=type)
new.shm.1 <- merge(new.shm, shm, by=c('clonotype.id', 'pos.nt', 'sample', 'type'))
df <- new.shm.1
save(df, file='sp_new.Rda')

