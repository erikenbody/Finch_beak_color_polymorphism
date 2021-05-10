library(data.table)
library(tidyverse)
library(zoo)
library(patchwork)
library(scales)

# whole genome ------------------------------------------------------------
chr_order3 <- c("chr1", "chr1A" , "chr2", "chr3", "chr4","chr4A", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11","chr12", "chr13", "chr14", "chr15", "chr17", "chr18", "chr19", "chr20", "chr21",
                "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chrZ",  "null")

l.scandens <- grep("ref.lrt.2.rec.lrt0.gz",grep("scandens_yellow_scandens_pink", list.files("data/gwas_angsd/Results_gwas_BALENCED/", recursive = T, full.names = T), value = T), value = T)

cat_scandens_lrt <- do.call("rbind",lapply(l.scandens,
                                      FUN=function(files){
                                        x <- fread(files, skip = 1)
                                        names(x) <- c("Chromosome","Position","Major","Minor", "Frequency","N","LRT","high_WT/HE/HO")
                                        x
                                      }))

cat_scandens_lrt$chr_ordered <- factor(cat_scandens_lrt$Chromosome, levels = chr_order3)
cat_scandens_lrt<-cat_scandens_lrt %>% arrange(chr_ordered, Position)
cat_scandens_lrt$row<-1:nrow(cat_scandens_lrt)

cat_scandens_lrt$fstrollmean <- zoo::rollmean(cat_scandens_lrt$LRT,100,fill=NA)
cat_scandens_lrt$chr_labels <- gsub("chr", "", cat_scandens_lrt$chr_ordered)

chr_breaks <- cat_scandens_lrt %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row))

png("output/gwas_angsd/scandens_manhattan_balanced.png", width = 16, height=1.75, res = 300, units = "in")
cat_scandens_lrt %>%  filter(LRT > -250) %>%
  #filter(LRT > 5) %>%
  ggplot(aes(x = row, y = LRT, col = chr_ordered)) + theme_bw()+
  scale_color_manual(values=rep(c("grey30","grey70"), length(levels(factor(cat_scandens_lrt$chr_ordered)))/2+1))+
  #scale_y_continuous(limits=c(0,0.12),breaks=seq(0,01,0.01), minor_breaks = NULL)+
  geom_point(size=1.5,shape=19,stroke=0.2) +
  scale_x_continuous(breaks = chr_breaks$chr_breaks, #https://stackoverflow.com/questions/50399838/how-to-alternate-a-new-line-for-overlapping-x-axis-labels
                     labels = function(labels) {
                       sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                     }) +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color = "black", vjust = 0.99, size = 10),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))
dev.off()


# fortis whole genome -----------------------------------------------------

l.fortis <- grep("ref.lrt.2.rec.lrt0.gz",grep("fortis_yellow_fortis_pink", list.files("data/gwas_angsd/Results_gwas_BALENCED/", recursive = T, full.names = T), value = T), value = T)

cat_fortis_lrt <- do.call("rbind",lapply(l.fortis,
                                           FUN=function(files){
                                             x <- fread(files, skip = 1)
                                             names(x) <- c("Chromosome","Position","Major","Minor", "Frequency","N","LRT","high_WT/HE/HO")
                                             x
                                           }))

cat_fortis_lrt$chr_ordered <- factor(cat_fortis_lrt$Chromosome, levels = chr_order3)
cat_fortis_lrt<-cat_fortis_lrt %>% arrange(chr_ordered, Position)
cat_fortis_lrt$row<-1:nrow(cat_fortis_lrt)

cat_fortis_lrt$fstrollmean <- zoo::rollmean(cat_fortis_lrt$LRT,100,fill=NA)
cat_fortis_lrt$chr_labels <- gsub("chr", "", cat_fortis_lrt$chr_ordered)

chr_breaks <- cat_fortis_lrt %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row))

png("output/gwas_angsd/fortis_manhattan_balanced.png", width = 16, height=1.75, res = 300, units = "in")
cat_fortis_lrt %>%  filter(LRT > -250) %>%
  #filter(LRT > 5) %>%
  ggplot(aes(x = row, y = LRT, col = chr_ordered)) + theme_bw()+
  scale_color_manual(values=rep(c("grey30","grey70"), length(levels(factor(cat_fortis_lrt$chr_ordered)))/2+1))+
  #scale_y_continuous(limits=c(0,0.12),breaks=seq(0,01,0.01), minor_breaks = NULL)+
  geom_point(size=1.5,shape=19,stroke=0.2) +
  scale_x_continuous(breaks = chr_breaks$chr_breaks, #https://stackoverflow.com/questions/50399838/how-to-alternate-a-new-line-for-overlapping-x-axis-labels
                     labels = function(labels) {
                       sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                     }) +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(color = "black", vjust = 0.99, size = 10),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))
dev.off()


# just chr24 ----------------------------------------------------

df.scandens.lrt2 <- fread("data/gwas_angsd/Results_gwas_BALENCED/scandens_yellow_scandens_pink %>% _BALENCED.ref.lrt.2.rec.lrt0.gz")
df.snp <- df.scandens.lrt2 %>% filter(Position == 6166878)

p.scandens.lrt2 <- df.scandens.lrt2 %>% filter(LRT!=-999.000000 ) %>%
  ggplot() + geom_point(aes(x = Position/1000000, y = LRT)) +
  geom_point(data = df.snp, aes(x = Position/1000000, y = LRT), color = "red") +
  theme(plot.title = element_text(face = "italic")) +
  theme_bw() + labs(title = "G. scandens", x = "Mb")
p.scandens.lrt2

df.scandens.lrt2.top10 <- df.scandens.lrt2 %>%  filter(LRT!=-999) %>% arrange(-LRT) %>% head(n = 10)

df.fortis.lrt2 <- fread("data/gwas_angsd/Results_gwas_BALENCED/fortis_yellow_fortis_pink_chr24_BALENCED.ref.lrt.2.rec.lrt0.gz")
df.snp <- df.fortis.lrt2 %>% filter(Position == 6166878)

p.fortis.lrt2 <- df.fortis.lrt2 %>% filter(LRT!=-999.000000 ) %>%
  ggplot() + geom_point(aes(x = Position/1000000, y = LRT)) +
  geom_point(data = df.snp, aes(x = Position/1000000, y = LRT), color = "red") +
  theme(plot.title = element_text(face = "italic")) +
  theme_bw() + labs(title = "G. fortis", x = "Mb")


df.fortis.lrt2.top10 <- df.fortis.lrt2 %>%  filter(LRT!=-999) %>% arrange(-LRT) %>% head(n = 10)

p.fortis.lrt2/p.scandens.lrt2
ggsave("output/gwas_angsd/fortis_scandens_LRT_recessive_chr24_balanced.png", width = 8, height = 6)

df.top10.fortis.scandens <- rbind(df.scandens.lrt2.top10, df.fortis.lrt2.top10)
df.top10.fortis.scandens %>% distinct(Position, .keep_all = T) %>%
  select(Chromosome, Position) %>% write_tsv("output/gwas_angsd/top_snps_fortis_scandens_balanced.txt", col_names = F)

tmp <- read.table("output/gwas_angsd/top_snps_fortis_scandens.txt")

# fortis and scandens combined  -------------------------------------------

df.combined.lrt2 <- fread("data/gwas_angsd/Results_gwas/Results_gwas/combined_yellow_combined_pink_chr24.res.lrt.2.rec.lrt0.gz")

df.combined.top100 <- df.combined.lrt2 %>% arrange(-LRT) %>%
  filter(Position < (6166878 + 200000) & Position > (6166878 - 200000)) %>%
  head(n = 100)
min(df.combined.top100$Position)
max(df.combined.top100$Position)
write_csv(df.combined.top100, "output/gwas_angsd/top100_snps_combined_gwas.txt")


p.combined.lrt2 <- df.combined.lrt2 %>% filter(LRT!=-999.000000 ) %>%
  ggplot() + geom_point(aes(x = Position/1000000, y = LRT), size = 3) +
  #geom_point(data = subset(df.combined.lrt2, Position %in% df.top10.fortis.combined$Position), aes(x = Position/1000000, y = LRT), color = "red") +
  theme(plot.title = element_text(face = "italic")) +
  theme_bw() +
  theme(text=element_text(size=30),
        panel.grid = element_blank(),
        plot.title = element_text(face = "italic")) +
  labs(title = "G. fortis & G. scandens", x = "chr24 (Mb)")

p.combined.lrt2
ggsave("output/gwas_angsd/fortis_scandens_LRT_recessive_chr24_COMBINED_balanced.png", width = 10, height = 10)

# mafs --------------------------------------------------------------------

df.fortis.yellow <- read.table(gzfile("data/gwas_angsd/Results_af/fortis_yellow_chr24_BALANCED.ref.mafs.gz"), header = T) %>%
  mutate(type = "fortis_yellow")
df.fortis.pink <- read.table(gzfile("data/gwas_angsd/Results_af/fortis_pink_chr24_BALANCED.ref.mafs.gz"), header = T) %>%
  mutate(type = "fortis_pink")
df.scandens.yellow <- read.table(gzfile("data/gwas_angsd/Results_af/scandens_yellow_chr24_BALANCED.ref.mafs.gz"), header = T) %>%
  mutate(type = "scandens_yellow")
df.scandens.pink <- read.table(gzfile("data/gwas_angsd/Results_af/scandens_pink_chr24_BALANCED.ref.mafs.gz"), header = T) %>%
  mutate(type = "scandens_pink")

df.all.af <- rbind(df.fortis.yellow, df.fortis.pink, df.scandens.yellow, df.scandens.pink)

df.all.af$type <- factor(df.all.af$type, levels = c("scandens_pink","fortis_pink", "scandens_yellow","fortis_yellow"))


df.all.af.wide <- df.all.af %>% select(position, unknownEM, type) %>%
  pivot_wider(id_cols = position, names_from = type, values_from = unknownEM) %>%
  mutate(mean_daf = abs(((fortis_yellow + scandens_yellow)/2) - ((fortis_pink + scandens_pink)/2)))

#remove any sites with missing data for one species and SNPs more than 400kb away
df.all.af.wide <- df.all.af.wide %>% filter(!is.na(fortis_yellow) & !is.na(fortis_pink) & !is.na(scandens_yellow) & !is.na(scandens_pink)) %>%
  filter(position < (6166878 + 400000) & position > (6166878 - 400000))

p.heatmap <- df.all.af %>%
  filter(position %in% df.all.af.wide$position) %>%
  #mutate(type = gsub("_"," ",type)) %>%
  ggplot() +
  geom_tile(aes(x = as.factor(position), y = type, fill = unknownEM)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low="pink", high="orange", mid = "white",
                       midpoint = (max(df.all.af$unknownEM, na.rm = TRUE) /2 )) +
  labs(fill = "AF", x = "Position", y = NULL)


p.daf <- ggplot(df.all.af.wide) + theme_bw() +
  geom_point(aes(x = factor(position), y = mean_daf)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "DAF")

p.daf / p.heatmap

layout <- "
A
B
B
B
"
p.daf + p.heatmap + plot_layout(design = layout)
ggsave("output/gwas_angsd/allele_frequency_top_SNPs.png", width = 10, height = 8)

df.exons <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_genes_exons.csv")
df.exons$gene_symbol <- toupper(df.exons$gene_symbol)

df.exons.overlap <- df.exons %>% filter(seqid == "24" & type == "exon" & start > min(df.all.af.wide$position) - 10000 & end < max(df.all.af.wide$position) + 10000)
df.exons.region <- df.exons %>% filter(gene_id %in% df.exons.overlap$gene_id & type == "exon") %>%
  mutate(gene_color = ifelse(gene_symbol == "BCO2","orange","not"))
df.genes.region <- df.exons %>% filter(gene_id %in% df.exons.region$gene_id & type == "gene") %>%
  mutate(gene_color = ifelse(gene_symbol == "BCO2","orange","not"))

p.daf2 <- ggplot() +
  geom_segment(data = df.genes.region, aes(y = 0.90, yend = 0.90, x = start, xend = end, color = gene_color), size = 1) +
  geom_segment(data = df.exons.region, aes(y = 0.90, yend = 0.90, x = start, xend = end, color = gene_color), size = 5) +
  geom_text(
    data = df.genes.region,
    mapping = aes(y = 0.95, x = ((start+end)/2), label = gene_symbol),
    parse = TRUE, hjust = 0, angle = 45, size = 3
  ) +
  geom_point(data = df.all.af.wide, aes(x = position, y = mean_daf)) + ylim(0, 1.05) + theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black","orange")) + ylab("DAF") + xlab(NULL) +
  scale_x_continuous(labels = label_number(scale = 1e-6)) + xlab("Mbp")

layout <- "
A
B
"
p.daf2 + p.heatmap + plot_layout(design = layout)

ggsave("output/gwas_angsd/allele_frequency_top_SNPs.pdf", width = 10, height = 8)
ggsave("output/gwas_angsd/allele_frequency_top_SNPs.png", width = 10, height = 8)
