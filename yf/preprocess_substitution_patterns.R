#Load data for both RNA and DNA

old <- c("Abdulain", "Ilgen", "Mamaev", "Smirnov", "Vlasov", "Xwoman")
young <- c("Antipyat", "Antipyat1", "Antipyat2", 
           "Epifancev", "Epifancev1", "Epifancev2",
           "Hadjibekov", "Koshkin", "Kovalchuk")

df <- data.frame()

load_sample <- function(sample, proj, type) {
  fileName <- paste('data/yf_', proj ,'_', type, "/", sample, ".shm.txt", sep = "")
  .df <- data.frame()
  
  if (file.exists(fileName)) {
    .df <- read.table(fileName, header=T, sep="\t")
    
    replica <- as.integer(substr(sample, nchar(sample), nchar(sample)))
    
    if (!is.na(replica)) {
      sample <- substr(sample, 1, nchar(sample) -1)
    } else {
      replica <- 1
    }
    
    .df$proj <- proj
    .df$type <- type
    .df$sample <- sample
    .df$replica <- replica
    
    .df <- subset(.df, mutation.type == "Substitution" & from.aa != "*" & to.aa != "*" & segment == "V")
  }
  
  .df
}

for (sample in old) {
  df <- rbind(df, load_sample(sample, "old", "DNA"))
  df <- rbind(df, load_sample(sample, "old", "RNA"))
}

for (sample in young) {
  df <- rbind(df, load_sample(sample, "young", "DNA"))
  df <- rbind(df, load_sample(sample, "young", "RNA"))
}

# Purge unused factor levels

df$segment.name <- factor(df$segment.name)
df$region <- factor(df$region)
df$from.nt <- factor(df$from.nt)
df$from.aa <- factor(df$from.aa)
df$to.nt <- factor(df$to.nt)
df$to.aa <- factor(df$to.aa)

# Summarize and remove alleles

library(plyr)
library(doMC)

doMC::registerDoMC(cores=50) 

df <- ddply(df, .(proj, sample, type, replica,
                  segment.name, region,
                  pos.nt, from.nt, to.nt,
                  pos.aa, from.aa, to.aa),
            summarize,
            rate = sum(freq),
            total.clonotypes = sum(clonotypes), 
            allele.rate = sum(freq / segment.freq) * length(segment.freq), 
            allele.clonotype.share = sum(clonotypes / segment.clonotypes) * length(segment.clonotypes),
            .parallel = T)

df <- subset(df, total.clonotypes <= 5 || allele.rate < 0.45 || allele.clonotype.share < 0.45)

save(df, file = "sp.Rda", compress = T)
