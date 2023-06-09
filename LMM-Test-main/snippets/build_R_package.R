
# build basic package skeleton in folder
devtools::create(".")

# documentation using roxygen2 package -> standard in R studio
# shortcut for building roxygen skeleton: Strf+Shift+Alt+R

# run documentation
devtools::document()

# install package 
devtools::install()

library(lme4randtest)

# for developers: load all functions into workspace
devtools::load_all()