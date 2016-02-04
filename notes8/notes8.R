## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----data_unadj----------------------------------------------------------
system("head unadjusted_unemployment.csv",intern=TRUE)
U1 <- read.table(file="unadjusted_unemployment.csv")
head(U1)

