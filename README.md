ResistanceGA
============

R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it it spossible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces.

To install this package, execute the following commands in R:

```
install.packages("devtools") # Installs the 'devtools' package
library(devtools) # Loads devtools

install_github("wpeterman/ResistanceGA") # Downloads package
require(ResistanceGA) # Installs package and the other required packages needed
```

###########################################

In order to use this package, you must have CIRCUITSCAPE installed.
Version 4.0 or higher is required:
Official releases can be found [here](https://code.google.com/p/circuitscape/downloads/list "CS downloads")


Further details of of how this method has been implemented can be found in:

Peterman, W. E., G. M. Connette, R. D. Semlitsch, and L. S. Eggert. in press. Ecological resistance surfaces predict fine scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology:10.1111/mec.12747.

Please cite this paper if you use this method!