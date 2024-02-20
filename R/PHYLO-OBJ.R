
# BIND DATA

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/MEIOFAUNA/INPUTS/"

# 1

tax <- list.files(path = wd, pattern = "CURATED-classify-consensus-blast-asvs.tsv", full.names = T)

ab_f <- list.files(path = wd, pattern = 'table_100_80', full.names = T)

ab <-  read_tsv(ab_f, skip = 1) %>% 
  dplyr::rename("Feature ID" = "#OTU ID") %>% 
  right_join(read_tsv(tax)) %>%
  # arrange("Feature ID") %>%
  data.frame(row.names = .$`Feature ID`) %>%
  select_if(is.double)

ab$Consensus <- NULL

dim(ab)

tax <-  read_tsv(tax) %>% arrange(match(`Feature ID`,rownames(ab) )) %>%  data.frame(row.names = .$`Feature ID`) 

tax$Feature.ID <- NULL

# 2
MTD <-  read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T)) %>% 
  select(`#SampleID`, Depth, Region) %>%
  mutate(Region = factor(Region, levels = c("Yucatan", "NW Shelf", "NW Slope", "Deep-sea"))) %>%
  dplyr::rename("LIBRARY_ID" = "#SampleID") %>%
  mutate(LIBRARY_ID = gsub("-", ".", LIBRARY_ID)) %>%
  filter(LIBRARY_ID %in% names(ab)) %>%
  data.frame(row.names = .$LIBRARY_ID)

identical(names(ab),rownames(MTD))

identical(rownames(ab), rownames(tax))

library(phyloseq)

phyloseq = phyloseq(otu_table(ab, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(MTD)
                    # phy_tree(tree)
                    )

saveRDS(phyloseq, paste0(wd, "phyloseq.rds"))
