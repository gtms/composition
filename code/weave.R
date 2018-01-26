library(rmarkdown)

render("../md/README.Rmd",
       ## output_dir = "../out/html",
       output_dir = "..",
       intermediates_dir = "../tmp")
