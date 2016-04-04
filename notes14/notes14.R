## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
  )

## ----prelims,echo=F,cache=F----------------------------------------------
set.seed(594709947L)
require(ggplot2)
theme_set(theme_bw())
require(plyr)
require(reshape2)
require(foreach)
require(doMC)
require(pomp)
stopifnot(packageVersion("pomp")>="0.69-1")

