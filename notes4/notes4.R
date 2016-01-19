## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----root----------------------------------------------------------------
roots <- polyroot(c(1,2,2))
roots

## ----abs_roots-----------------------------------------------------------
abs(roots)

