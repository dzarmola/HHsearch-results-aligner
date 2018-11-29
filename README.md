# HHsearch-results-aligner
Creates a plot visualizing multiple sequence alignment based on pairwise HHsearch results

HHsearch results (.hhr files) can either be specified as positional arguments, or taken from a directory (--hhr_dir - this will glob this directory with "/\*.hhr").

Cluster file is an optional MCL dump file which will add additional lines to the plot which indicate cluster allegiance
for each profile.

Default e-value (above which HHsearch alignment will be disregarded) is 0.001, it can be change with --eval

Example plot:
![Example plot](relative/path/to/img.jpg?raw=true "Alignment with clustering info")


From help:
usage: my_little_merger_error_cor.py [-h] [--hhr_dir HHR_DIR] [--eval EVAL]
                                     [--plot_from_save PLOT_FROM_SAVE]
                                     [--cluster_file CLUSTER_FILE]
                                     [--save_name SAVE_NAME]
                                     [--plot_name PLOT_NAME]
                                     [HHR [HHR ...]]

Merge hhpred alignments

positional arguments:
  HHR                   .hhr files to be merged

optional arguments:
  -h, --help            show this help message and exit
  --hhr_dir HHR_DIR     Directory from which all .hhr will be taken
  --eval EVAL           Maximum e-value of used hhsearch results
  --plot_from_save PLOT_FROM_SAVE
                        Savepoint from a previous run
  --cluster_file CLUSTER_FILE
                        Cluster defining file: each cluster in one line
  --save_name SAVE_NAME
                        Custom savepoint name
  --plot_name PLOT_NAME
                        Custom plot name
