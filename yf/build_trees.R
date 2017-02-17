library(stringr)
library(dplyr)
library(fitdistrplus)
library(doParallel)
library(foreach)
library(igraph)

find_kmers <- function(string, k=5){
  n <- nchar(string) - k + 1
  kmers <- str_sub(string, 1:n, 1:n + k - 1)
  return(kmers)
}

info <- function(kmer){
  return(kmer.df$neg_ln[kmer.df$kmer == kmer])
}

sum_info <- function(l){
  return( sum(unlist(lapply(l, info))) )
}

mutations_weight <- function(mutations){
  # unused function in current version
  w = 1
  #mutations = str_split(str, ',')[[1]]
  for (m in mutations){
    mut = str_extract(m, '[>TAGC]+')
    if (str_sub(m, 1, 1)  == 'S'){
      w = w * filter(sub.df, from.to==mut)$freq
    }
    else{
      w = w * filter(indel.df, type==str_sub(m, 1, 1), len==length(mut))$freq
    }
  }
  return(w)
}

inter_info <- function(i){
  inf = data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("n1", "n2", "mutual", "inf1", "inf2", "mutations"))),
                  stringsAsFactors=F)
  for (j in (i+1):length(kmers)){
    if (df$v[i] == df$v[j] & df$j[i] == df$j[j]){
      mut.i = df$all.mutations[[i]]
      mut.j = df$all.mutations[[j]]
      shared_mutations = intersect(mut.i, mut.j)
      
      if (length(shared_mutations) == length(mut.i) | length(shared_mutations) == length(mut.j)){
        n1 = ifelse(length(shared_mutations) == length(mut.i), i, j)
        n2 = ifelse(n1 == i, j, i)
        shared_kmers = intersect(kmers[[i]], kmers[[j]])
        mutations_between = setdiff(df$all.mutations[[n2]], df$all.mutations[[n1]])
        
        row <- list(n1 = n1, n2 = n2, mutual = sum_info(shared_kmers), 
                    inf1 = sum_info(kmers[[n1]]), inf2 = sum_info(kmers[[n2]]),
                    mutations = NA, shared_kmers_num = length(shared_kmers),
                    w = mutations_weight(mutations_between), mut_between = length(mutations_between))
        inf <- rbind(inf, row)
        inf$mutations[nrow(inf)] <- paste(mutations_between, collapse=',')
      }
    }
  }
  return(inf)
}

interleave <- function(v1, v2){
  z <- c()
  for (i in 1:length(v1)){
    z <- c(z, v1[i], v2[i])
  }
  return(z)
}

tree_statistics <- function(tree){
  clone <- data.frame(ndn = character(), freq = double(), freq_sum = double(), leaves = integer(), 
                      nodes = integer(), mut_in_root = integer(), diameter = integer(), 
                      mean_mut = double(),total_mut = integer(), branching = double())
  if (length(V(tree)) > 1){
    edges <- get.edgelist(tree)
    root <- setdiff(edges[,1], edges[,2])
    leaves <- setdiff(edges[,2], edges[,1])
    freq <- df[as.integer(root),]$freq
    ndn <- df[as.integer(root),]$ndn
    mut_in_root <- df$all.mutations[[as.numeric(root)]]
    v <- df[as.integer(root),]$v
    j <- df[as.integer(root),]$j
    
    mut_from_root <- c()
    for (l in leaves){
      mut_in_leave = df$all.mutations[[as.numeric(l)]]
      shared_mutations = intersect(mut_in_root, mut_in_leave)
      mut_from_root <- c(mut_from_root, length(mut_in_leave) - length(shared_mutations))
    }
    diameter = max(mut_from_root)
    mean_mut = mean(mut_from_root)
    mut_sum = sum(mut_from_root)
    
    freq_sum = sum(df[as.numeric(names(V(tree))),]$freq)
    path_length <- sapply(shortest_paths(tree, root)$vpath, length)
    leaves_n = length(leaves)
    nodes_n = length(V(tree))
    branching = leaves_n/mean(path_length)
    clone <- rbind(clone, list(ndn = ndn, freq = freq, v=v, j=j, freq_sum = freq_sum, leaves = leaves_n, nodes = nodes_n, 
                mut_in_root = length(mut_in_root), diameter = diameter, mean_mut = mean_mut,
                total_mut = mut_sum, branching = branching))
  }
  return(clone)
}




old <- paste0('yf_old_RNA/', c("Abdulain", "Ilgen", "Mamaev", "Smirnov", "Vlasov"))
young <- paste0('yf_young_RNA/', c("Antipyat", "Epifancev", "Hadjibekov", "Koshkin", "Kovalchuk"))

for (sample in c(old, young)){
  #prepare dataset
  df <- read.table(paste0('~/yf/', sample, '.clones.txt'), header=T, sep="\t")
  df <- df %>% mutate(all.mutations=paste(mutations.nt.FR1, mutations.nt.CDR1, mutations.nt.FR2, mutations.nt.CDR2,
                           mutations.nt.FR3, mutations.nt.CDR3, mutations.nt.FR4, sep=','),
                      cdr3nt = as.character(cdr3nt)) %>%
    mutate(all.mutations = str_extract_all(all.mutations, '([\\w\\d>:]+)'),
           ndn = str_sub(cdr3nt, pmax(0, v.end.in.cdr3-4), pmin(j.start.in.cdr3+4, nchar(cdr3nt))))
  #df <- df %>% group_by(cdr3nt, v.end.in.cdr3, j.start.in.cdr3, v, j, all.mutations)  %>% summarise(sample = 'raji')
  df$v <- str_sub(str_extract(df$v, '(.+)\\*'), 1, -2)
  df$j <- str_sub(str_extract(df$j, '(.+)\\*'), 1, -2)
  
  #import frequencies of k-mers and mutations and fit_gamma
  load('~/yf/frequencies.rda')
  kmer.df$neg_ln <- -log(kmer.df$freq)
  sub.df$neg_ln <- -log(sub.df$freq)
  indel.df$neg_ln <- -log(indel.df$freq)
  
  # list all 5-mers
  cl <- makeCluster(8)
  registerDoParallel(cl)
  kmers <- foreach(i = df$ndn, .packages='stringr') %dopar% find_kmers(i)
  
  # get and filter edges
  inter <- foreach(x = 1:(length(kmers)-1), .combine='rbind', .packages = c('dplyr', 'stringr')) %dopar% inter_info(x)
  
  inter <- inter %>% mutate(p_value = 1 - pgamma(mutual/(inf1*inf2), shape = fit_gamma$estimate['shape'], rate = fit_gamma$estimate['rate']),
                            edge = paste(n1,n2,sep='_')) %>% filter(p_value < 1e-08)
  
  final_edges <- foreach(x = 1:length(kmers), .combine='rbind', .packages = c('dplyr')) %dopar% get_arbor_edges(x)
  
  #write net file
  net <- data.frame(from = final_edges$n1, to = final_edges$n2, interaction = rep('shm', nrow(final_edges)))
  write.table(net, file=paste0('~/yf/trees/', sample, '.net.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  # write edge file
  edge <- data.frame(edge = paste(final_edges$n1, '(shm)', final_edges$n2, sep=' '), mutations = final_edges$mutations)
  write.table(edge, file=paste0('~/yf/trees/', sample, '.edge.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  # write node file
  all_nodes <- unique(c(final_edges$n1, final_edges$n2))
  node <- df[all_nodes,] %>% dplyr::select(cdr3aa, v, j, freq)
  node$node <- all_nodes
  write.table(node, file=paste0('~/yf/trees/', sample, '.node.txt'), sep='\t', row.names=FALSE, quote=FALSE)
  
  #graph analysis
  g <- graph( edges=interleave(final_edges$n1, final_edges$n2) )
  g <- set.vertex.attribute(g, 'name', value = as.character(1:length(V(g))))
  g <- set.edge.attribute(g, 'weight', value = final_edges$mut_between)
  components <- decompose.graph(g)
  clones <- foreach(x = components, .combine='rbind', .packages = c('igraph')) %dopar% tree_statistics(x)
  singletons <- df[setdiff(1:nrow(df), all_nodes), ] %>% dplyr::select(ndn, v, j, freq, all.mutations) %>%
    mutate(freq_sum = freq, leaves = 0, nodes = 1, diameter = 0,
           mean_mut = 0, total_mut = 0, branching = 0)
  mut_in_root <- sapply(singletons$all.mutations, length)
  singletons <- dplyr::select(cbind(singletons, mut_in_root), -all.mutations)
  clones <- rbind(clones, singletons)
  write.table(clones, file=paste0('~/yf/trees/stat/', sample, '.txt'), sep='\t', row.names=FALSE, quote=FALSE)
}
