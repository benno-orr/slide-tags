# Source another Rmd's code chunks as plain R (purl -> source).
rmd = function(rmd, rmd_dir="/n/data1/hms/scrb/chen/lab/bco/scripts/rmd") {
    tf <- tempfile(fileext = ".R")
    knitr::purl(file.path(rmd_dir, paste0(rmd, ".Rmd")), output = tf, documentation = 0, quiet = TRUE)
    source(tf, local = FALSE)
    print(paste0("loaded ", rmd))
}
