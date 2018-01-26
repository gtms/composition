library(rmarkdown)

render("../md/readme.Rmd",
       ## output_dir = "../out/html",
       output_dir = "..",
       intermediates_dir = "../tmp")
