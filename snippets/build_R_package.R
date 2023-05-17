
# build basic package skeleton in folder
devtools::create(".")

# documentation using roxygen2 package -> standard in R studio

# run documentation
devtools::document()

# install package 
devtools::install()

library(lme4randtest)
?