library(devtools)
wd <- "/Users/priscilla/Documents/Documents - Priscilla’s MacBook Pro/repositories/MixedPsy"
setwd(wd)
document()
setwd("..")
install("MixedPsy")
setwd(wd)
library(MixedPsy)
