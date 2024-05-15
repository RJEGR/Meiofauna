# After upgrade DIVERSITY.R
# Correlate
# DIVERSITY2COR.R


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)
library(phyloseq)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"


R <- sample_data(read_rds(paste0(wd, "/ps.rds"))) %>% as(., "matrix")

vars_to_numeric <- c("Depth","Latitude",	"Longitude", 
  "Clay",	"Silt",	"Sand",	"IC",	"TOC",	"CN",	"Oxygen", 
  "Finas",	"Medias",	"Muy_finas",	"Gruesas")


measures <- c("Observed", "Chao1", "Shannon", "InvSimpson", "Evenness")

reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")

R <- R %>% as_tibble() %>% mutate_at(c(vars_to_numeric, measures), as.numeric) 


R %>%
  pivot_longer(cols = measures, values_to = "Diversity") %>%
  filter(name == "Shannon") %>%
  mutate(Region = factor(Region, levels = reg_levels)) %>%
  mutate(name = factor(name, levels = measures))%>% 
  # mutate(x = TOC) %>%
  ggplot(aes(x = TOC, y = Diversity)) +
  facet_wrap(~ name, scales = "free_y", nrow = 1) +
  geom_smooth(method = "lm", se = T, color = 'black', formula = y ~ x) +
  geom_point() +
  ggpubr::stat_cor(color = 'black')

# ggpmisc::stat_poly_eq(aes(label = stat(eq.label)), formula = y ~ x, parse = T)

# multiple corr-mat


library(rstatix)

# R %>%
#   drop_na(Oxygen) %>%
#   rstatix::cor_mat(vars = c("Shannon", "Oxygen", "TOC", "Latitude", "Longitude")) %>%
#   cor_reorder() %>%
#   pull_lower_triangle() %>%
#   cor_plot(label = TRUE)


measures <- names(R)[names(R) %in% measures] 


cor_out <- vector("list", length(vars_to_numeric))

names(cor_out) <- vars_to_numeric

for (i in 1:length(vars_to_numeric)) {
  
  var <- vars_to_numeric[i] # vars_to_numeric[1]
  
  which_vars <- c(measures, var)
  
  cor_out[[i]] <- R %>% mutate(param = var) %>% drop_na(param) %>% rstatix::cor_mat(vars = which_vars) %>% cor_gather()
  
  
}

do.call(bind_rows, cor_out) -> cor_df 

cor_df %>% count(var1)
cor_df %>% count(var2)


cor_df %>% filter(var1 %in% measures & var2 %in% measures & cor != 1)

cor_df <- cor_df %>% filter(!var1 %in% measures & cor != 1) 

cor_df <- cor_df %>%
  mutate(star = ifelse(p <.001, "***", 
    ifelse(p <.01, "**",
      ifelse(p <.05, "*", ""))))

lo = floor(min(cor_df$cor))
up = ceiling(max(cor_df$cor))
mid = (lo + up)/2

cor_df %>%
  mutate(var1 = factor(var1, levels = vars_to_numeric)) %>%
  # mutate(cor = ifelse(p <.05, NA, cor)) %>%
  ggplot(aes(y = var1, x = var2, fill = cor)) +
  geom_tile(color = 'black', linewidth = 0.7, width = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  geom_text(aes(label = star), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Richness compared to multiple c. subjects") +
  theme_minimal(base_size = 14, base_family = "GillSans") +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
      margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm"))) -> p


ggsave(p, filename = 'alpha_diversity_cor.png', path = wd, width = 3.5, height = 5, device = png, dpi = 300)


cor_df %>%
  filter(var1 == "Shannon") %>% 
  filter(var2 != "Shannon") %>% 
  distinct(cor, .keep_all = T) %>%
  # mutate(y_star = cor + (0.2+SE)*sign(logFC))
  arrange(desc(cor), .by_group = T) %>% 
  mutate(Label = var2, row_number = row_number(Label)) %>%
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  ggplot(aes(x = cor, y = Label)) +
  # geom_col(fill = "white", color = "black", linewidth = 1.5) 
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.5) +
  geom_text(aes(label=star), vjust = -0.25, color="black", 
    position=position_dodge(width = .5), size = 4, family = "GillSans") +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  labs(x = "Correlation", y = "Shannon compared to multiple subjects", title = "") +
  theme_minimal(base_size = 14, base_family = "GillSans") +
  theme(plot.title = element_text(hjust = 0.5))
  

# R %>%
#   pivot_longer(cols = c(measures, "Evenness"), values_to = "rich_val", names_to = "rich_name") %>%
#   filter(rich_name %in% "Shannon") %>%
#   pivot_longer(cols = c(measures, "Evenness"), values_to = "rich_val", names_to = "rich_name")
#   # mutate(Region = factor(Region, levels = reg_levels)) %>%
#   # mutate(name = factor(name, levels = c(measures, "Evenness"))) %>%
#   group_by(Batch) %>% cor_test(vars = "Nominal", vars2 = c("Input"))


# scatterplot of signif


filt_var <- cor_df %>% filter(star != "") %>% distinct(var1) %>% pull()

R %>%
  select(c(measures, vars_to_numeric)) %>%
  pivot_longer(cols = vars_to_numeric, names_to = "vars", values_to = "val_var") %>%
  drop_na(val_var) %>%
  filter(vars %in% filt_var) %>%
  # pivot_longer(cols = vars_to_numeric, names_to = "vars") %>%
  # mutate(name = factor(name, levels = measures))%>% 
  # mutate(x = TOC) %>%
  ggplot(aes(x = val_var, y = Shannon)) +
  facet_wrap(~ vars, scales = "free_x", nrow = 1) +
  geom_smooth(method = "lm", se = T, color = 'black', formula = y ~ x) +
  geom_point(alpha = 0.5, size = 0.5) +
  ggpubr::stat_cor(color = 'black')
