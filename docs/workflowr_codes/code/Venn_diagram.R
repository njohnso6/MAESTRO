install.packages("VennDiagram")
library(VennDiagram)
library(tidyverse)
library(tidyverse)
frag <- read_tsv("./data/fragments/fragments_corrected_count.tsv", col_names = FALSE)
frag
frag <- frag %>% mutate(coord = str_c(X1,X2,X3, sep = "_"))
set1 <- frag$coord
####################

cellranger_frag <- read_tsv("cellranger_frag.tsv", col_names = FALSE)
minimap_frag <- read_tsv("minimap_frag.tsv", col_names = FALSE)
chromap_frag <- read_tsv("chromap_frag.tsv", col_names = FALSE)

cellranger_frag <- cellranger_frag %>% mutate(coord = str_c(X1,X2,X3, sep = "_"))
minimap_frag <- minimap_frag %>% mutate(coord = str_c(X1,X2,X3, sep = "_"))
chromap_frag <- chromap_frag %>% mutate(coord = str_c(X1,X2,X3, sep = "_"))

set1 <- cellranger_frag$coord
set2 <- minimap_frag$coord
set3 <- chromap_frag$coord
###################

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("cellranger" , "minimap2" , "chromap"),
  filename = 'frag_venn_diagramm.png',
  output=TRUE
)




# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
myCol <- brewer.pal(3, "Accent")
myCol_3 <- brewer.pal(3, "Set3")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("cellranger (10,889,126)" , "minimap2 (7,593,511)" , "chromap (10,754,168)"),
  filename = 'frag_venn_diagramm_5.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 600 ,
  width = 600 ,
  resolution = 500,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,

  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("cellranger (10,889,126)" , "minimap2 (7,593,511)" , "chromap (10,754,168)"),
  filename = 'frag_venn_diagramm_5.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  fill = c(alpha("#0073C2FF",0.3), alpha('#EFC000FF',0.3), alpha('#CD534CFF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  rotation = 1)
