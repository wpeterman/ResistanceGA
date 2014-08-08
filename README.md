ResistanceGA
============

### An R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it is possible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces. Resistance distances can be calculated as cost distances (least cost path) between points, or as circuit-based resistance distances calculated using CIRCUITSCAPE    

To install this package, execute the following commands in R:

```
# Install 'devtools' package, if needed
if(!("devtools" %in% list.files(.libPaths()))) {
    install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE) 
} 

library(devtools) # Loads devtools

install_github("wpeterman/ResistanceGA") # Download package
library(ResistanceGA) # Installs package and the other required packages needed
```
Once the package is installed, you can further explore the functions by opening the HTML 'Vignette' using the code below, or you can view the vignette directly [**here**](https://dl.dropboxusercontent.com/u/23513016/ResistanceGA.html "Vignette")
```
vignette('ResistanceGA')  # Opens tutorial in web browser
```
*****

### Other notes

If you wish to optimize using CIRCUITSCAPE, you must have CIRCUITSCAPE installed.
Version 4.0 or higher is required.
Official CIRCUITSCAPE releases can be found [**here**](https://code.google.com/p/circuitscape/downloads/list "CS downloads")

At this point in time, CIRCUITSCAPE can only be executed through `ResistanceGA` on **_Windows_** machines. The only modifications that would be necessary to execute these functions on another platform would be a determining how to execute CIRCUITSCAPE from R. If there is interest in being able to use these functions on other platforms, I can work to incorporate that capability in future releases.


This approach has been developed from the methods first utilized in Peterman et al. (2014). I have also written a preprint manuscript further describing this `R` package. **_Please cite this paper if you use these methods!_**

Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402–2413. [**PDF**](http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/peterman_et_al._2014--mec.pdf "Peterman et al.")

Peterman, W.E. 2014. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. bioRxiv doi: 10.1101/007575. [**PDF**](http://biorxiv.org/content/early/2014/07/29/007575 "Peterman bioRxiv")