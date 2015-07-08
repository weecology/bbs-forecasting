library(prism)
options(prism.path = "./data/prismdata")
get_prism_monthlys(type="tmean", year = 1990:2010, month = 6, keepZip=F)

datapath = "./data/prismdata/"
datadirs = dir(datapath)
for (datadir in datadirs) {
  bil_file = paste(datadir, '.bil', sep = "")
  bil_file_path = file.path(datapath, datadir, bil_file)
  sql_file = paste(datadir, ".sql", sep = "")
  sql_file_path = file.path(datapath, datadir, sql_file)
  system(paste("raster2pgsql -s 4326", bil_file_path, datadir, ">", sql_file_path))
  system(paste("psql -d bbsforecasting -f", sql_file_path))
}