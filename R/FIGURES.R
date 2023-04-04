# DATAVIZ

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# 0) LOAD DATA =====

wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"

subdir <- paste0(wd, "CURATED_DB_DIR")

tax_f <- read_tsv(list.files(path = subdir, pattern = 'taxonomy.tsv$', full.names = T))

ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)

MTD <-  read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T))

hist(tax_f$Confidence)

# 1) SPLIT TAX ======

into <- paste0("rank", 0:9)

tax_f <- tax_f %>%
  select(-`Feature ID` , -Confidence) %>%
  separate(Taxon, sep = ";", into = into) # %>% 
  # mutate_at(into, list(~ str_replace_all(., c("os__" = ""))))

tax_f %>% view()
