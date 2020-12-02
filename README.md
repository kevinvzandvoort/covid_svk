# IBM model code for "The effectiveness of population-wide, rapid antigen test based screening in reducing SARS-CoV-2 infection prevalence in Slovakia"

## Dependencies

### R packages:
```R
  install.packages("Rcpp", "RcppGSL", "data.table", "ggplot2", "patchwork", "readxl", "distcrete", "qs", "tictoc")
```

### C++ libraries
  ```cpp
  #include <iostream>
  #include <vector>
  #include <cstdlib>
  #include <Rcpp.h>
  #include <omp.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>
```

* To run model scenarios, use script: `scenario_modelling.R`
** calls `./scripts/functions.R`
** calls `./scripts/read_data.R`
** compiles `./model/covidIBM.cpp`

* To create plots from output, use script: `process_outputs.R`
