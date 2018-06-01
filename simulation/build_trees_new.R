library(stringr)
library(dplyr)
library(doParallel)
library(foreach)
library(igraph)
library(muscle)
library(Biostrings)
library(tidyr)
library(purrr)
library(data.table)
library(stringi)

select = dplyr::select
summarise = dplyr::summarise

skf <- function(x, K){
  #resulting columns are ordered as k-mers in kmer.df
  s <- stri_count_fixed(x, pattern = filter(kmer.df, k==K)$kmer, overlap=TRUE)
  as.data.frame(t(s))
}

count_kmers <- function(df, K, full){
  if (full){
    kmers <- do.call(rbind, mclapply(df$ndn, function(x) skf(x, K), mc.cores=50))
  } else {
    kmers <- do.call(rbind, mclapply(df$cdr3nt, function(x) skf(x, K)))
  }
  colnames(kmers) <- filter(kmer.df, k==K)$kmer
  kmers
}

JSD <- function(c.table, weight, pseudocount = 1e-6){
  m <- as.data.frame(t(c.table))
  m[m==0] <- pseudocount
  colnames(m) <- c('c1', 'c2')
  m$w = weight
  
  m = m %>%
    mutate(f1 = c1/sum(c1), 
           f2 = c2/sum(c2), 
           mean = 0.5 * (f1 + f2),
           kld1 = f1 * log2(f1 / mean) / w,
           kld2 = f2 * log2(f2 / mean) / w)
  
  0.5 * (sum(m$kld1) + sum(m$kld2))
}

m_dist <- function(c.table, bg){
  k1 = colnames(c.table)[c.table[1,] != 0]
  k2 = colnames(c.table)[c.table[2,] != 0]
  shared.kmers = colnames(c.table[,colSums(c.table != 0) == 2])
  
  if (length(shared.kmers) == 0){
    return(NA)
  } else {
    shared.kmer.inf = sum(unlist(sapply(shared.kmers, function(x) -log(bg$freq[bg$kmer == x]))))
    inf1 = sum(unlist(sapply(k1, function(x) -log(bg$freq[bg$kmer == x]))))
    inf2 = sum(unlist(sapply(k2, function(x) -log(bg$freq[bg$kmer == x]))))
    m = shared.kmer.inf/(inf1+inf2)
    log(m) 
  }
}

kmer_dist <- function(i, table, method, bg.freq){
  df = data.frame()
  for (j in (i+1):nrow(table)){
    print(j)
    if (method == 'JSD'){
      dist = JSD(table[c(i,j),], weight = rep(1, nrow(bg.freq)))
    } else if (method == 'JSDW'){
      dist = JSD(table[c(i,j),], weight = bg.freq$freq)
    } else if (method == 'm'){
      dist = m_dist(table[c(i,j),], bg = bg.freq)
    }
    df <- rbind(df, data.frame(n1 = i, n2 = j, dist = dist))
  }
  df
}

find_root <- function(tree, dt){
  all_nodes <- V(tree)$name
  nodes <- dt[all_nodes]
  nodes[which.min(mut.num)]$clone
}

set_edge_direction <- function(tree, root){
  bfs <- bfs(tree, root=root, order=TRUE, father=TRUE)
  f = as.vector(bfs$father) # igraph indices of vertices' fathers, not names
  is.not.root = !is.na(f) #remove NA produced by tree root
  n1 = get.vertex.attribute(tree, 'name',  f[is.not.root])
  n2 = V(tree)$name[is.not.root]
  data.table(n1 = n1, n2 = n2)
}

mutations_weight <- function(mutations){
  prod(merge(data.frame(mutation = mutations), frequencies)$log.freq)
}

get_mutations <- function(s1, s2){
  DNAset = DNAStringSet(c(s1, s2))
  alignment <- unmasked(muscle(DNAset, quiet = TRUE))
  aln.table <- align2table(alignment) %>%
    mutate(s1 = as.character(s1), s2 = as.character(s2))
  
  ind <- which(aln.table$s1 != aln.table$s2)
  m <- aln.table[ind,] %>%
    filter(s1 != '-', s2 != '-') %>%
    mutate(c1 = paste(s1,s2,sep='>'), c2 = paste(s2,s1,sep='>')) %>%
    select(-s1, -s2)
  m
}

get_parent <- function(m1, m2){
  score1 <- mutations_weight(m1)
  score2 <- mutations_weight(m2)
  parent = ifelse(score1 > score2, 'c1', 'c2')
  parent
}



##TEST

setwd('simulation')
#import frequencies of k-mers and mutations
load('mut_frequencies.rda')
load('kmer_freq.Rda')

opt <- list(a = 0.45, f = T)

load('mixcr_1000_6_0.25.Rda')
#df <- filter(df, root < 20)

df <- df %>%
  mutate(contignt = as.character(contignt),
         ndn = str_sub(contignt, vEnd+1, jStart-1)) %>%
  filter(nchar(ndn) > 4)
df$all.mutations = sapply(df$all.mutations, function(x) str_split(x, ',')[[1]])
df$mut.num = sapply(df$all.mutations, length)
df$mut.idx = sapply(df$all.mutations, function(x) str_match(x, '(.+):')[,2])
df$clone = as.character(df$clone) #to avoid malignant behaviour of igraph and data.table keys
dt <- as.data.table(df)

#cl <- makeCluster(detectCores())
cl <- makeCluster(50)
registerDoParallel(cl)


y <- df %>%
  filter(parent %in% clone) %>%
  filter(parent!=clone) %>%
  select(clone, parent) %>%
  mutate(pair = paste(parent, clone, sep='_')) %>%
  arrange(parent, clone)

vj <- unique(select(df, v, j))

for (t in c(-6)){
  opt$a <- t
  edges <- data.table()

  for (x in 1:nrow(vj)){
    cat(t, x, "VJ pair out of", nrow(vj), "(", round(x/nrow(vj)*100, 2), "%)\n")
    cat(vj[x,]$v, vj[x,]$j, '\n')
    vt = vj[x,]$v
    jt = vj[x,]$j
    .dt.vj <- dt[v == vt & j == jt]
    setkey(.dt.vj, clone)
    print(nrow(.dt.vj))
    
    kmers <- count_kmers(.dt.vj, K=4, full=T)
    bg.freq <- filter(kmer.df, k==4) #add here VJ specific kmer frequencies?
    
    # for (i in 1:nrow(.dt.vj)){
    #   print(i)
    #   kmer_dist(i, kmers, method = 'JSD', bg.freq)
    # }
    
    dist <- foreach(x=1:(nrow(kmers)-1), 
                    .combine='rbind',
                    .packages = c('dplyr')) %dopar%
      kmer_dist(x, kmers, method = 'm', bg.freq)
    
    #return global clone numbers instead of .dt.vj numbers
    dist$n1 <- .dt.vj[dist$n1]$clone
    dist$n2 <- .dt.vj[dist$n2]$clone
  
    pairs <- dist %>% filter(dist > opt$a)
  
    mut_num_between <- function(x, mode='v_m'){
      #v_i takes reverse mutations into account (though prediction quality is the same)
      #v_m is 2 times faster
      if (mode == 'v_m'){
        m1 <- .dt.vj[pairs$n1[x]]$all.mutations[[1]]
        m2 <- .dt.vj[pairs$n2[x]]$all.mutations[[1]]
        length(c(setdiff(m1, m2), setdiff(m2, m1)))
      } else if(mode == 'v_i'){
        m1 <- .dt.vj[pairs$n1[x]]$all.mutations[[1]]
        m2 <- .dt.vj[pairs$n2[x]]$all.mutations[[1]]
        i1 <- .dt.vj[pairs$n1[x]]$mut.idx[[1]]
        i2 <- .dt.vj[pairs$n2[x]]$mut.idx[[1]]
        i1 <- i1[which(m1 %in% setdiff(m1,m2))]
        i2 <- i2[which(m2 %in% setdiff(m2,m1))]
        length(unique(c(i1, i2)))
      }
    }
    
    .edges <- data.frame()
    if (nrow(pairs)>0){
      pairs$mut.num <- sapply(1:nrow(pairs), function(x) mut_num_between(x, mode='v_m'))
      
      pairs <- pairs %>% mutate(n1 = as.character(n1), n2 = as.character(n2))
      g <- graph_from_edgelist(as.matrix(pairs %>% select(n1, n2)), directed = F)
      #E(g)$weight <- pairs$dist
      E(g)$weight <- pairs$mut.num
      components <- decompose.graph(g)
      
      .edges <- foreach(x = components, .combine='rbind', .packages=c('igraph', 'data.table')) %dopar% {
        tree <- mst(x) # get minimal spanning tree
        root <- find_root(tree, .dt.vj) #find root as vertex closest to germline
        set_edge_direction(tree, root) #root the tree
      }
      
      edges <- rbind(edges, .edges)
    }

    cat(nrow(.edges), 'new edges,', nrow(edges), 'total edges\n')
  }
  
  fit <- edges %>%
    mutate(pair = paste(n1, n2, sep='_')) %>%
    arrange(n1, n2)
  
  metrics(y$pair, fit$pair)
}


