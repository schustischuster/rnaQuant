# rnaQuant
RNA quantification from tissue sections.

## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and visualization](#data-analysis-and-visualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install and load the following R packages before running the reproducible script:

```R
# Required packages
lib_List <- c("dplyr", "ggplot2", "grid", "scales")

# Install missing packages
instpack <- lib_List %in% installed.packages()[,"Package"]
if (any(instpack == FALSE)) {
  install.packages(lib_List[!instpack])
}

# Load packages
invisible(lapply(lib_List, library, character.only = TRUE))

```

### Data input
Download and extract the rnaQuant repository to the working directory on your computer. Then, set the path for input and output files and source the R scripts:

```R
in_dir <- file.path("rnaQuant-master", "data")
out_dir <- file.path("rnaQuant-master")
path_to_R_files <- file.path("rnaQuant-master", "R")

# Source R files
sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
   }
}
 
sourceDir(path_to_R_files)
```
---
## Data analysis and visualization

...Add data description and processing here...


* `plotRIN.R(temp, ...)`


| Arguments  |  |
| :---  | :---  |
| temp  | Indicates data set. Values can be either "56" or "62". |


To reproduce the results of this study, execute the following function calls:

```R
plotRIN(temp = "62")
plotRIN(temp = "56")

```

## Session info
