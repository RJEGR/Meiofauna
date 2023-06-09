#!/usr/bin/env Rscript

# El siguiente script identifica el primer campo (columna, como el identificador de la especie [processid]); 
# así como la última posición (ie. última columna) a la secuencia de origen [voucher_type];
# Si existen identificadores repetidos, el script re-etiquetara con el sufijo equivalente a las n repeticiones del identificador separado por un guión bajo “_”. 
# El resto de las columnas son implementadas para formatear la taxonomía para el algoritmo de mothur (classify.seqs). 
# El nombre del archivo (Ej. taxonomy.Animals.csv) es utilizado para etiquetar el segundo nivel taxonómico (Ej. root;rank1, donde rank1 pertenece al nivel taxonómico reino) por lo que es necesario nombrar el archivo de entrada con la etiqueta correspondiente. 
# Los campos de la asignación taxonómica vacíos son rellenados con valores NA, y finalmente se retienen sólo los taxones con resolución hasta especie. 

# How tu run:
# bold_Plants_Chlorophyta.trim, bold_Plants_Magnoliophyta.trim, bold_Plants_Bryophyta.trim
# ~$ kingdom='Plants'
# ~$ cut -f1,3,5,7,9,13,15,17 *${kingdom}*trim | sed 's/ /_/g' > taxonomy.${kingdom}.csv 
# Rscript --vanilla bold_public_process_for_RDP.R taxonomy.${kingdom}.csv
args = commandArgs(trailingOnly=TRUE)

cat('1. Reading file at time:', (d <- format(Sys.time(), "%H:%M:%S")), '\n')


# file = args[1]

# cat  SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.tax | awk '{gsub(/;/, "\t" , $2); print $0}'  > SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.tax2
# EL METODO DE AWK NO ESTA LIMPIANDO TODAS LAS LECTURAS. PROBAR QIIMEE PROX SEMANA DESDE CERO

file = list.files(path = "/Users/cigom/Documents/MEIOFAUNA_PAPER/mothur/", 
  pattern = "SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.tax2", full.names = T)

# By default keep species is TRUE

keep_spp <- TRUE

kingdom <- strsplit(file, '[.]')[[1]][2]
root.rank <- paste0('root;Eukaryota',';', kingdom)
colnames <- c('Id', 'Domain', 'Clade1','Clade2','Phyla','Clade3','Phyla2','subkingdom','Clase', 'Orden', 'Familia', 'Genero', 'Especie')
# read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
x <- read.csv(file, header=FALSE, sep='\t', na.strings=c("","NA"), stringsAsFactors=FALSE)

library(tidyverse)

x %>% separate(V1, sep = " ", into = c("ID", "Domain")) -> x

sort(table(x$Domain))

# fungus                      Archaea 
# 191                          192                         1699 
# Eukaryota                     Bacteria 
# 20744                        97326 

which_dom <- c("fungus", "Archaea", "Eukaryota", "Bacteria")

x %>% filter(Domain %in% which_dom) %>% count(Domain)
x %>% filter(Domain %in% which_dom) %>% count(V2, sort = T)

names(x) <- colnames

cat('2. Keeping only the Specie resolution ...\n')


if (keep_spp) { x <- x[complete.cases(x),] }




Id <- make.unique(as.vector(x[,1]), sep = "_")

cat('3. Adding root and ...\n')

y <- cbind(root.rank, (x[-c(1,ncol(x))]))


cat('4. Formating taxonomy ...\n')

library(tidyr)

save <- cbind(Id, unite(y, sep = ";", remove = TRUE, col = 'Taxonomy'))

save$Taxonomy <- sapply(save$Taxonomy, 
  function(x){gsub(pattern = "$",
    replacement = ";", x)})


cat('5. Formating fasta ...\n')

n_seqs <- length(Id)
seq_headers <- vector(n_seqs, mode="character")

for (i in 1:n_seqs) {
  seq_headers[i] <- paste(">", Id[i], sep = "")
}

fasta <- c(rbind(seq_headers, x[,ncol(x)]))

cat('6. Writing outputs ... \n')
# 1
write.table(save, file = paste0("Bold.",kingdom, ".tax"), append = FALSE, quote = FALSE, sep = " ",
  na = "NA", row.names = FALSE,
  col.names = FALSE)
# 2
write(fasta,
  file=paste0("Bold.",
    kingdom, 
    ".fasta"))

cat('7. DONE!\n')

difftime(strptime(format(Sys.time(), "%H:%M:%S"), format = "%H:%M:%S"), 
  strptime(d, format = "%H:%M:%S"),units="auto")

quit(save ='no')

