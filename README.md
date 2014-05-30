ResistanceGA
============

<<<<<<< HEAD
An R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it is possible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces.
=======
R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it is possible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces.
>>>>>>> 3d5c85b0f6c4899222388e6be9fa35b662b3f685

To install this package, execute the following commands in R:

```
install.packages("devtools") # Installs the 'devtools' package
library(devtools) # Loads devtools

install_github("wpeterman/ResistanceGA") # Downloads package
require(ResistanceGA) # Installs package and the other required packages needed
```
<<<<<<< HEAD
Once the package is installed, you can further explore the functions by opening the HTML 'Vignette'. You can also view the vignette directly [*here*](https://dl.dropboxusercontent.com/u/23513016/ResistanceGA_Vignette.html "Vignette")
```
vignette('ResistanceGA_Vignette')  # Opens tutorial in web browser
```
=======
Once the package is installed, you can further explore the functions by opening the HTML 'Vignette' locatedin the help files.
>>>>>>> 3d5c85b0f6c4899222388e6be9fa35b662b3f685

### Other notes

In order to use this package, you must have CIRCUITSCAPE installed.
Version 4.0 or higher is required:
Official CIRCUITSCAPE releases can be found [here](https://code.google.com/p/circuitscape/downloads/list "CS downloads")
<<<<<<< HEAD
=======

>>>>>>> 3d5c85b0f6c4899222388e6be9fa35b662b3f685


This approach has been developed from the methods first utilized in Peterman et al. (2014). Please cite this paper if you use these methods!

Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402–2413. [PDF](http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/peterman_et_al._2014--mec.pdf "Peterman et al.").