N_files = 500
N_jobs = 25

files_per_job = N_files / N_jobs

starts = seq(0, N_jobs - 1) * files_per_job + 1
ends = starts + files_per_job - 1

script = "#!/bin/bash
  
# Job name and who to send updates to
#SBATCH --job-name=mistnet
#SBATCH --mail-user=harris.d@ufl.edu
#SBATCH --mail-type=FAIL,END

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
#SBATCH --mem-per-cpu=8G

# Job run time in [DAYS]
# HOURS:MINUTES:SECONDS
# [DAYS] are optional, use when it is convenient
#SBATCH --time=30-00:00:00

# Save some useful information to the 'output' file
date;hostname;pwd

# Load R and run a script
module load R
Rscript --default-packages=stats,graphics,grDevices,utils,methods analysis/fit-mistnet.R"

for (i in 1:N_jobs) {
  filename = paste0("job_", i, ".job")
  # send the script to SLURM with the specified starts & ends as arguments
  cat(script, starts[[i]], ends[[i]], "\n", file = filename)
  system(paste("sbatch", filename), wait = FALSE)
  print(paste("job", i, "started"))
  file.remove(filename)
}
