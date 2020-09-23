ResistanceGA
============

### An R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it is possible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces. Resistance distances can be calculated as cost distances (least cost path) between points, or as circuit-based resistance distances calculated using CIRCUITSCAPE    

To install this package, execute the following commands in R:

```
# Install 'devtools' package, if needed
if(!("devtools" %in% list.files(.libPaths()))) {
    install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE) 
} 

devtools::install_github("wpeterman/ResistanceGA", build_vignettes = TRUE) # Download package

library(ResistanceGA) # Installs package and the other required packages needed
```
Once the package is installed, you can further explore the functions by opening the 'Vignette' using the code below, or you can view the vignette directly [**here**](http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/resistancega.pdf "Vignette")
**NOTE: The Vignette is quite dated, and has not been updated to fully reflect the features now incorporated in the current version.**
```
vignette('ResistanceGA')  # Opens PDF
```
*****

### Other notes

Optimization with CIRCUITSCAPE (v4) is still possible, although not actively supported in the most recent version.
If you wish to optimize using CIRCUITSCAPE, you it is highly recommended that you install Julia and the CIRCUITSCAPE Julia package
General instructions [**here**](https://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/julia_guide.pdf "Julia Guide")


***
This approach has been developed from the methods first utilized in Peterman et al. (2014). I have also written a preprint manuscript further describing this `R` package. **_Please cite these papers if you use these methods!_**

Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402–2413. [**PDF**](http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/peterman_et_al._2014--mec.pdf "Peterman et al.")

Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. Methods in Ecology and Evolution doi:10.1111/2041-210X.12984. [**PDF**](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12984 "MEE Publication")
