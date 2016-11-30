library(plyr)
library(ggplot2)
library(dplyr)
library(reshape2)

old_dna <- c("Abdulain", "Ilgen", "Mamaev", "Smirnov", "Vlasov", "Xwoman")
young_dna <- c("Antipyat1", "Antipyat2", "Epifancev1", "Epifancev2", "Hadjibekov", "Koshkin", "Kovalchuk")
old_rna <- c("Abdulain", "Ilgen", "Mamaev", "Smirnov", "Vlasov")
young_rna <- c("Antipyat", "Epifancev", "Hadjibekov", "Koshkin", "Kovalchuk")
path = 'immunity/analysis(2.11)/data/'

df <- data.frame()

for (sample in old_rna) {
  .df <- read.table(paste(path, 'yf_old_RNA/', sample, ".clones.txt", sep = ""), header=T, sep="\t")
  .df$proj <- "old"
  .df$sample <- sample
  df <- rbind(df, .df)
}

for (sample in young_rna) {
  .df <- read.table(paste(path, 'yf_young_RNA/', sample, ".clones.txt", sep = ""), header=T, sep="\t")
  .df$proj <- "young"
  .df$sample <- sample
  df <- rbind(df, .df)
}

for (sample in 1:9) {
  .df <- read.table(paste('immunity/analysis(2.11)/data/normal/age_ig_s', sample, '_RO.txt', sep = ""), header=T, sep="\t")
  .df$proj <- "normal"
  .df$sample <- sample
  df <- rbind(df, .df)
}

df.1.v <- ddply(df, .(v, proj, sample, cdr3aa), summarize, n=1)
df.2.v <- ddply(df.1.v, .(v, proj, sample), summarize, n=sum(n))
df.2.v <- ddply(df.2.v, .(proj, sample), transform, share = n/sum(n))
df.3.v <- dcast(df.2.v, sample + proj ~ v, value.var = 'share')
df.3.v[is.na(df.3.v)] <- 0
rownames(df.3.v) <- df.3.v$sample
df.3.v <- subset(df.3.v, select = -c(sample, proj))
proj_col <- ifelse(rownames(df.3.v) %in% old_rna, 'red', ifelse(rownames(df.3.v) %in% young_rna, 'blue', 'black'))
heatmap.2(as.matrix(df.3.v), RowSideColors = proj_col, scale="row")

#Skip low frequensy genes
x <- df.3.v
for (v in names(x)){
  if (sum(x[v]) < 0.07){
    x[v] <- NULL
  }
}
heatmap.2(as.matrix(x), RowSideColors = proj_col, scale="row")

ggplot(df.3.v, aes(x = v, y = share, fill = proj)) + geom_boxplot() + theme_bw() + xlab('')

#N insertion length
df.1 <- mutate(df, ndn = j.start.in.cdr3 - v.end.in.cdr3)
df.1 <- mutate(df.1, nn = (d.start.in.cdr3 - v.end.in.cdr3) + (j.start.in.cdr3-d.end.in.cdr3))
ggplot(df.1[!(is.na(df.1$d)),], aes(x = nn, color=proj)) + geom_density()
ggplot(df.1, aes(x = ndn, color=proj)) + geom_density()

ggplot(df.1,aes(x=nn)) + 
  geom_histogram(data=subset(df.1,proj == 'old'),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(df.1,proj == 'young'),fill = "blue", alpha = 0.2

ggplot(df.1,aes(x=ndn)) + 
  geom_histogram(data=subset(df.1,proj == 'old'),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(df.1,proj == 'young'),fill = "blue", alpha = 0.2)