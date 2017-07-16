CV = FALSE

if (CV) {
  N_files = 250
} else {
  N_files = 500
}
N_jobs = 250

files_per_job = N_files / N_jobs

stopifnot(files_per_job == as.integer(files_per_job))

starts = seq(0, N_jobs - 1) * files_per_job + 1
ends = starts + files_per_job - 1

# The "script" object has the whole job script EXCEPT the name of the .R file
# to run and the arguments that get sent to it.
script = "#!/bin/bash
  
# Job name and who to send updates to
#SBATCH --job-name=mistnet
#SBATCH --mail-user=harris.d@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

# Where to put the outputs:
#   %A expands to the job-name specified above
#   %j expands into the job number (a unique identifier for this job)
#SBATCH --output mistnet-%j.out
#SBATCH --error mistnet-%j.err

# Number of nodes to use, following https://wiki.rc.ufl.edu/doc/R#FAQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Memory per cpu core. Default is megabytes, but units can be specified 
# with M or G for megabytes or Gigabytes.
#SBATCH --mem-per-cpu=6G

# Job run time in [DAYS]
# HOURS:MINUTES:SECONDS
# [DAYS] are optional, use when it is convenient
#SBATCH --time=96:00:00

# Save some useful information to the 'output' file
date;hostname;pwd

# Load R and run a script
module load R
Rscript --default-packages=stats,graphics,grDevices,utils,methods"

for (i in 1:N_jobs) {
  filename = ifelse(CV, "analysis/cv-mistnet.R", "analysis/fit-mistnet.R")
  jobname = paste0("job_", i, ".job")
  # send the script to SLURM with the specified filename, starts & ends as 
  # arguments, plus N_files (used for reproducibility in the CV script)
  cat(script, filename, starts[[i]], ends[[i]], N = N_files, "\n", file = jobname)
  system(paste("sbatch", jobname), wait = TRUE)
  file.remove(jobname)
  Sys.sleep(0.25)
}
