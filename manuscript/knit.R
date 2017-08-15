library(knitr)
library(rmarkdown)
library(tidyverse)

# Eliminate newlines inside of R code marked by `r code_goes_here`.
readLines("manuscript/draft.Rmd") %>% 
  paste(collapse = "\n") %>% 
  gsub("`r\n", "`r ", .) %>% 
  cat(file = "manuscript/draft-modified.Rmd")

render("manuscript/draft-modified.Rmd", 
       output_format = "pdf_document", 
       output_dir = "manuscript",
       output_file = "draft.pdf")
file.remove("manuscript/draft-modified.Rmd")
