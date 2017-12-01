library(rmarkdown)

render("../md/supplementary.Rmd",
       output_dir = "../out/html",
       intermediates_dir = "../tmp")
