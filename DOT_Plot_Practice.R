#Trying to plot my paf file 

paf <- read.table(
  "magister_vs_floridiana.paf",
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = "",
  comment.char = ""
)

plot(paf$V8, paf$V3,
     pch = ".",
     cex = 0.5,
     col = rgb(0,0,0,0.2),
     xlab = "Target position",
     ylab = "Query position",
     main = "Raw PAF dot plot")

paf <- read.table("magister_vs_floridiana.paf",
                  sep="\t",
                  header=FALSE,
                  fill=TRUE,
                  quote="")

query_start  <- paf$V3
target_start <- paf$V8

plot(target_start, query_start,
     pch=".", cex=0.5)
