library(stringr)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)

find_kmers <- function(string, k=5){
  n <- nchar(string) - k + 1
  kmers <- str_sub(string, 1:n, 1:n + k - 1)
  return(kmers)
}

info <- function(kmer){
  return(-log(kdf$neg_lnp[kdf$kmer == kmer]))
}

sum_info <- function(l){
  return( sum(unlist(lapply(l, info))) )
}

#prepare datasets

df <- read.table('migec.txt', header=T, sep="\t")
df <- mutate(df, ndn = str_sub(CDR3.nucleotide.sequence, Last.V.nucleotide.position-4, First.J.nucleotide.position+4))
#df <- read.table('yf/age_ig/age_ig_s1_RO.txt', header=T, sep="\t")
#df <- mutate(df, ndn = str_sub(cdr3nt, v.end.in.cdr3-3, j.start.in.cdr3+3))
df <- df[sample(nrow(df)),]
df$sample = rep('migec', nrow(df))

df2 <- read.table('raji_R12.txt', header=T, sep="\t")
df2 <- df2[sample(nrow(df2)),]
df2 <- df2 %>% group_by(cdr3nt, v.end.in.cdr3, j.start.in.cdr3, v, j)  %>% summarise(sample = 'raji') %>%
  mutate(ndn = str_sub(cdr3nt, v.end.in.cdr3-4, j.start.in.cdr3+4))
df2$v <- str_sub(str_extract(df2$v, '(.+)\\*'), 1, -2) #because mutate does not work that's why
df2$j <- str_sub(str_extract(df2$j, '(.+)\\*'), 1, -2) 

df3 <- read.table('age_ig/age_ig_s1_RO.txt', header=T, sep="\t")
df3 <- rbind(df3, read.table('age_ig/age_ig_s2_RO.txt', header=T, sep="\t"))
df3 <- rbind(df3, read.table('age_ig/age_ig_s3_RO.txt', header=T, sep="\t"))
df3 <- df3 %>% group_by(cdr3nt, v.end.in.cdr3, j.start.in.cdr3, v, j)  %>% summarise(sample = 'age')
df3 <- mutate(df3, ndn = str_sub(cdr3nt, v.end.in.cdr3-4, j.start.in.cdr3+4))
df3$v <- str_sub(str_extract(df3$v, '(.+)\\*'), 1, -2) #because mutate does not work that's why
df3$j <- str_sub(str_extract(df3$j, '(.+)\\*'), 1, -2) 
df3 <- df3[sample(nrow(df3)),]
  
df4 <- data.frame(ndn = c(df$ndn[1:3000], df2$ndn, df3$ndn[1:3000]), 
                  sample = c(df$sample[1:3000], df2$sample, df3$sample[1:3000]),
                  v = c(as.character(df$V.segments)[1:3000], df2$v, df3$v[1:3000]),
                  j = c(as.character(df$J.segments)[1:3000], df2$j, df3$j[1:3000]))
df4$ndn <- as.character(df4$ndn)

# list all 5-mers
cl <- makeCluster(8)
registerDoParallel(cl)
kmers <- foreach(i = df4$ndn, .packages='stringr') %dopar% find_kmers(i)

# count k-mers
all_kmers <- unlist(kmers)
kdf <- data.frame(kmer = all_kmers)
kdf <- kdf %>% group_by(kmer) %>% summarise(count = n())
kdf <- mutate(kdf, neg_lnp = -log(count/sum(count)))
hist(kdf$count)

inter_info <- function(i){
  inf <- data.frame(n1 = integer(0), n2 = integer(0), sample1 = character(0), sample2 = character(0), 
                    mutual = double(0), inf1 = double(0), inf2 = double(0))
  for (j in (i+1):length(kmers)){
    if (df4$v[i] == df4$v[j] & df4$j[i] == df4$j[j]){
      row <- list(n1 = i, n2 = j, sample1 = df4$sample[i], sample2 = df4$sample[j], 
                  mutual = sum_info(intersect(kmers[[i]], kmers[[j]])), 
                  inf1 = sum_info(kmers[[i]]), inf2 = sum_info(kmers[[j]]))
      inf <- rbind(inf, row)
    }
  }
  return(inf)
}

inter <- foreach(x = 1:(length(kmers)-1), .combine='rbind', .packages = 'stringr') %dopar% inter_info(x)
stopCluster(cl)

inter = mutate(inter, samples = ifelse(sample1 == 'raji' & sample2 == 'raji', 'raji',
                                       ifelse(sample1 == 'migec' & sample2 == 'migec', 'migec', 
                                              ifelse(sample1=='age' & sample2 == 'age', 'age', 'diff'))))

ggplot(inter[!is.na(inter$samples),]) + geom_density(aes(x=mutual, color = samples))
ggplot(inter[!is.na(inter$samples),]) + geom_histogram(aes(x=mutual, color = samples), binwidth=1, position='dodge') + scale_y_log10()
