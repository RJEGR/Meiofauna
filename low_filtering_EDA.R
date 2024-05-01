
# EDA
# Filtering out low abundance asvs (spurious)
# Prepare data
# Since we are interested in alpha diversity, it is probably not a bad idea to prune OTUs that are not present in any of the samples 

tax <- as(tax_table(phyloseq), "matrix") %>% as_tibble(rownames = "Feature ID")

Freq <- function(x){rowSums(x > 1)}

prevelancedf <- as(otu_table(phyloseq), "matrix") %>%
  as_tibble(rownames = "Feature ID") %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.integer))),
    Prevalence = Freq(across(where(is.integer)))) %>%
  select(`Feature ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) %>%
  left_join(tax)


# out %>% view()

n_sam <- ncol(as(otu_table(phyloseq), "matrix"))

prevelancedf %>% 
  # drop_na(k) %>%
  ggplot(aes(TotalAbundance, Prevalence/n_sam)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  # scale_color_manual(values = c('red', 'blue')) +
  # labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Consensus", caption = caption_1) +
  facet_wrap( ~ p) +
  theme_bw(base_size = 17, base_family = "GillSans") +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> p2

p2 + theme(panel.border = element_blank(), legend.position = "top")

# Now lets investigate low prevelance/abundance phylum and subset them out.

summary_prevalence <- plyr::ddply(prevelancedf, "p", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

# summary_prevalence %>% arrange(desc(total_abundance)) %>% view()

# Using the table above, determine the phyla to filter based on the 0.001 threshold

sum(summary_prevalence$total_abundance)*0.001
table(summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= 0.001)

keepPhyla <- summary_prevalence$p[summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= 0.001]

physeq = subset_taxa(phyloseq, p %in% keepPhyla)

summary_prevalence <- summary_prevalence[summary_prevalence$p %in% keepPhyla,]

summary_prevalence

# Individual taxa filtering
# Subset to the remaining phyla by prevelance.

prevelancedf1 = subset(prevelancedf, p %in% get_taxa_unique(physeq, taxonomic.rank = "p"))

ggplot(prevelancedf1, aes(TotalAbundance,Prevalence / nsamples(physeq),color=p)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.10, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~p) + theme(legend.position="none")

# here we can see a break above the 0.25 % with exception in proteobacteria

#  Define prevalence threshold as 10% of total samples ~ set of replicates
prevalenceThreshold = 0.10 * nsamples(physeq)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevelancedf1)[(prevelancedf1$Prevalence >= prevalenceThreshold)]

length(keepTaxa)

physeq2 <- prune_taxa(keepTaxa, physeq)
physeq2

# Agglomerate taxa at the Genus level (combine all with the same name) keeping all taxa without genus level assignment

length(get_taxa_unique(physeq, taxonomic.rank = "p"))
physeq2_glom = tax_glom(physeq, "p", NArm = FALSE)

physeq2_glom

sum(colSums(otu_table(physeq2_glom))) /sum(colSums(otu_table(phyloseq)))

# Now lets filter out samples (outliers and low performing samples)
# Do some simple ordination looking for outlier samples, first we variance stabilize the data with a log transform, the perform PCoA using bray’s distances



logt  = transform_sample_counts(physeq2_glom, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "MDS", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples",
  color = "Tissue", shape = "Time") + 
  labs(col = "Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))
