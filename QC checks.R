library(GWASinspector)

## import the QC-configuration file 
job <- setup_inspector("GWASinspector/config.ini")

## check the created instance
## input result files that will be inspected are also displayed
job

## run the algorithm 
job <- run_inspector(job)

## check the results
## comprehensive report and result file are already saved in the output folder
result_inspector(job)
