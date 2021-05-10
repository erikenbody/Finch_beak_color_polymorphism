
df.allele.share <- read.table('data/allele_sharing/sum_alt_all_species.txt')

#df.allele.share <- df.allele.share %>% filter(V1!=0)
head(df.allele.share$V1)

df.sum.share <- df.allele.share %>% group_by(V1) %>% summarise(n.sp = n()) 
  

ggplot() + geom_col(data = df.sum.share, aes(x = V1, y = n.sp)) + 
  geom_text(data = df.sum.share, aes(x = V1, y = n.sp + 150000, label = n.sp), size = 3) +
  theme_bw() + xlab("Number of species") + ylab("Number of shared SNPs") 
#dir.create("output/allele_sharing")
ggsave("output/allele_sharing/plot_allele_sharing.png", height = 6, width = 8)
ggsave("output/allele_sharing/plot_allele_sharing.pdf", height = 6, width = 8)

#percent shared all sp
df.sum.share[15,2] / sum(df.sum.share$n.sp) *100
#percent shared all but 1 sp
df.sum.share[14,2] / sum(df.sum.share$n.sp) *100

#14 + 15 species
df.sum.share[15,2] / sum(df.sum.share$n.sp) *100 + df.sum.share[14,2] / sum(df.sum.share$n.sp) *100
