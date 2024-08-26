library(phyloseq)
library(tidyverse)
library(vegan)
library(DESeq2)
library(ANCOMBC)
library(ComplexHeatmap)


#Build phyloseq project
otu <- read.table(file = "feature-table.tsv", sep = "\t", header = T, row.names = 1,
                  skip = 1, comment.char = "")
taxonomy <- read.table(file = "taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# clean the taxonomy, Greengenes format
tax <- taxonomy %>%
  select(Taxon) %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

metadata <- read.table(file = "sample-metadata.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
TREE = read_tree("tree.nwk")
# merge the data
ps <- phyloseq(OTU, TAX, SAMPLE,TREE)

#Check the sequencing depth with rarefaction curves

set.seed(111)
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=1103, replace=F)

#Alpha diversity

png(filename="Alpha_Diversity_Plot.png")
plot_richness(ps.rarefied, x="body.site", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
dev.off()

rich = estimate_richness(ps.rarefied, measures = c("Observed", "Shannon"))
wilcox.observed <- pairwise.wilcox.test(rich$Observed,
                                        sample_data(ps.rarefied)$body.site,
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

wilcox.shannon <- pairwise.wilcox.test(rich$Shannon,
                                       sample_data(ps.rarefied)$body.site,
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

#Beta diversity
#PCoA 
dist = phyloseq::distance(ps.rarefied, method="bray")
ordination = ordinate(ps.rarefied, method="PCoA", distance=dist)
png(filename="Beta_Diversity_Plot_1.png")
plot_ordination(ps.rarefied, ordination, color="body.site") +
  theme_classic() +
  theme(strip.background = element_blank())
dev.off()

#PERMANOVA/ADONIS
metadata <- data.frame(sample_data(ps.rarefied))
test.adonis <- adonis(dist ~ body.site, data = metadata)
test.adonis <- as.data.frame(test.adonis$aov.tab)

#PAIRWISE PERMANOVA
cbn <- combn(x=unique(metadata$body.site), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
metadata_sub <- data.frame(sample_data(ps.subs))
permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ body.site,
                             data = metadata_sub)
p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)

#NMDS and ANOSIM
dist = phyloseq::distance(ps.rarefied, method="jaccard", binary = TRUE)
ordination = ordinate(ps.rarefied, method="NMDS", distance=dist)
png(filename="Beta_Diversity_Plot_2.png")
plot_ordination(ps.rarefied, ordination, color="body.site") +
  theme_classic() +
  theme(strip.background = element_blank())
dev.off()

#ANOSIM
metadata <- data.frame(sample_data(ps.rarefied))
anosim(dist, metadata$body.site)

#PAIRWISE ANOSIM
cbn <- combn(x=unique(metadata$body.site), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
metadata_sub <- data.frame(sample_data(ps.subs))
permanova_pairwise <- anosim(phyloseq::distance(ps.subs, method="jaccard", binary = TRUE),
                             metadata_sub$body.site)
p <- c(p, permanova_pairwise$signif[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)

#Abundance bar plot
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(body.site, Phylum) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,body.site,Phylum) %>%
  summarise(Abundance=sum(Abundance))

png(filename="Abundance_Bar_Plot_Phylum.png")
ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", aes(fill=Phylum)) +
  labs(x="", y="%") +
  facet_wrap(~body.site, scales= "free_x", nrow=1) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = -90))
dev.off()

ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Genus', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Genus <- as.character(ps.melt$Genus)

ps.melt <- ps.melt %>%
  group_by(body.site, Genus) %>%
  mutate(median=median(Abundance))
# select group mean > 1
keep <- unique(ps.melt$Genus[ps.melt$median > 2.5])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 2.5%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,body.site,Genus) %>%
  summarise(Abundance=sum(Abundance))

png(filename="Abundance_Bar_Plot_Genus.png")
ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus)) + 
  labs(x="", y="%") +
  facet_wrap(~body.site, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))
dev.off()

#Differential abundance analysis
sample_data(ps)$body.site <- as.factor(sample_data(ps)$body.site) # factorize for DESeq2
ps.taxa <- tax_glom(ps, taxrank = 'Species', NArm = FALSE)

# pairwise comparison between gut and tongue
ps.taxa.sub <- subset_samples(ps.taxa, body.site %in% c("gut", "tongue"))
# filter sparse features, with > 90% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)
ps_ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ body.site)
# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(ps_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
# Only keep gut and tongue samples
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.pse.sub)), ps.taxa.rel.sig)

matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
# Define the annotation color for columns and rows
annotation_col = data.frame(
    Subject = as.factor(metadata_sub$subject),
    `Body site` = as.factor(metadata_sub$body.site),
    check.names = FALSE
)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
    Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
)
rownames(annotation_row) = rownames(matrix)

# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
    Subject = c(`subject-1` = "red", `subject-2` = "blue"),
    `Body site` = c(gut = "purple", tongue = "yellow"),
    Phylum = phylum_col
)

png(filename="HeatMap_DESeq2.png")
ComplexHeatmap::pheatmap(matrix, scale= "row",
                         annotation_col = annotation_col,
                         annotation_row = annotation_row,
                         annotation_colors = ann_colors)
dev.off()

#ANCOM-BC
# pairwise comparison between gut and tongue
ps.taxa.sub <- subset_samples(ps.taxa, body.site %in% c("gut", "tongue"))
# ancombc
out <- ancombc(phyloseq = ps.taxa.sub, formula = "body.site", 
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
              group = "body.site", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res <- out$res
# select the bottom 20 with lowest p values
res.or_p <- rownames(res$q_val["body.sitetongue"])[base::order(res$q_val[,"body.sitetongue"])]
taxa_sig <- res.or_p[1:20]
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
# Only keep gut and tongue samples 
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.sub)), ps.taxa.rel.sig)

matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
# Define the annotation color for columns and rows
annotation_col = data.frame(
    Subject = as.factor(metadata_sub$subject),
    `Body site` = as.factor(metadata_sub$body.site),
    check.names = FALSE
)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
    Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
)
rownames(annotation_row) = rownames(matrix)

# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
    Subject = c(`subject-1` = "red", `subject-2` = "blue"),
    `Body site` = c(`gut` = "purple", tongue = "yellow"),
    Phylum = phylum_col
)
png(filename="HeatMap_ANCOM.png")
ComplexHeatmap::pheatmap(matrix,
                         annotation_col = annotation_col,
                         annotation_row = annotation_row,
                         annotation_colors = ann_colors)
dev.off()


 






