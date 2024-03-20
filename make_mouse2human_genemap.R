### This script converts HOM_MouseHumanSequence into a dictionary table of mouse-human gene symbols 
# deleted row 2714 (duplicate)
library(tidyverse)
source('code-dependencies/qtl_functions.R')

db <- read_csv('source_data/HOM_MouseHumanSequence.csv')
names(db)[1:2] <- c('DB_Class_Key', 'Common_Organism_Name')
db <- db[,c('DB_Class_Key', 'Common_Organism_Name', 'Symbol')]

# mapdb <- pivot_wider(db, id_cols = DB_Class_Key, names_from = Common_Organism_Name, values_from = Symbol)
# names(mapdb) <- c('DB_Class_Key', 'mouse', 'human')
# mapdb <- unnest(mapdb, cols = c('mouse', 'human'))

mapdb <- pivot_wider(db, names_from = Common_Organism_Name, values_from = Symbol, values_fn = function(x) paste(x, collapse=","))
names(mapdb) <- c('DB_Class_Key', 'mouse', 'human')

# separate longer so that each row has only one mouse gene and one human gene
mapdb <- separate_longer_delim(mapdb, cols = 'human', delim = ',')

ensure_directory('derived_data')
write_csv(mapdb, 'derived_data/mouse_human_gene_map.csv')
