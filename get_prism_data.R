library(prism)
options(prism.path = "./data/prismdata")
get_prism_monthlys(type="tmean", year = 1990:2010, month = 6, keepZip=F)
