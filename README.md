# ResistanceGA

### An R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it is possible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces. Resistance distances can be calculated as cost distances (least cost path) between points, or as circuit-based resistance distances calculated using CIRCUITSCAPE

To install this package, execute the following commands in R:

    # Install 'devtools' package, if needed
    if(!("devtools" %in% list.files(.libPaths()))) {
        install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE) 
    } 

    devtools::install_github("wpeterman/ResistanceGA", build_vignettes = TRUE) # Download package

    library(ResistanceGA) # Loads package and the other dependnecies

If installation with `build_vignettes = TRUE` fails, try installing the `tinytex` R package and then attempt to install `ResistanceGA` again.

Once the package is installed, you can view the 'Getting Started' vignette in R.

------------------------------------------------------------------------

### Other notes

Optimization with CIRCUITSCAPE (v4) is still possible, although not actively supported in the most recent version. If you wish to optimize using CIRCUITSCAPE, it is highly recommended that you install Julia and the CIRCUITSCAPE Julia package. General instructions [**here**](https://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/julia_guide.pdf "Julia Guide"). There is also a `Julia_Guide` vignette with the package now.

------------------------------------------------------------------------

This approach has been developed from the methods first utilized in Peterman et al. (2014). The first formal analysis using ResistanceGA was Ruiz-López et al. (2016). The primary citation for the package is Peterman (2018), Methods in Ecology.

Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402–2413. [**PDF**](http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/peterman_et_al._2014--mec.pdf "Peterman et al.")

Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. Methods in Ecology and Evolution 9, 1638–1647.
<doi:10.1111/2041-210X.12984>. [**PDF**](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12984 "MEE Publication")

Ruiz-López, M.J., Barelli, C., Rovero, F., Hodges, K., Roos, C., Peterman, W.E., Ting, N., 2016. A novel landscape genetic approach demonstrates the effects of human disturbance on the Udzungwa red colobus monkey (*Procolobus gordonorum*). Heredity, Ruiz-López 116, 167–176. <https://doi.org/10.1038/hdy.2015.82>
