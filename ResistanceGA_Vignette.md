Example Function Use: ResistanceGA
=====
## An R Package for Optimizing Resistance Surfaces using Genetic Algorithms

With this vignette/tutorial, hopefully you'll get an idea of what each of the functions in this package can do, as well present an example (using simulated data) of how one can optimize resistance surfaces in isolation as well as simultaneously to create novel resistance surfaces. This 'package' (I use that term very loosely) has largely been developed from functions I wrote to conduct different landscape genetic analyses. See [Peterman et al. (2014)](http://onlinelibrary.wiley.com/doi/10.1111/mec.12747/abstract "Published Molecular Ecology Study") for the original conception of optimizing resistance surfaces using optimization functions. This approach was limited to optimization of continuous surfaces in isolation. Since that paper, I've further developed the optimization method to utilize gentic algorithms, as implements using the `ga` function from the [GA package](http://cran.r-project.org/web/packages/GA/index.html "GA package, CRAN") in R. By moving to genetic algorithms, much more comlex parameter space can be effectively searched, which allowed for the optimization of categorical resistance surfaces, as well as optimization of multiple resistance surfaces simultaneously.

A few words of caution. I have made every effort to run and test each function with simulated data, but I make no guarantees concerning function performance and stability. Data formatting can be a challenge, and I have tried to simplify the process as much as possible. If errors occur, start by making sure that you are providing function inputs in the correct format. Hopefully this tutorial will help with that process. Lastly, this is not a fast process. Even with the 80x80 pixel simulated landscapes used in this tutorial, each optimization iteration takes 0.75--1.25 seconds to complete (Intel i7 3.4 GHz processor, 24 BG RAM). The largest surfaces I've attempted to optimize using these methods was 600x600 pixels, which took ~13 seconds per iteration. Depending upon whether you are optimizing a single surface or multiple surfaces simultaneously, the genetic algorithms typically run for 50--150 generations. `ga` settings will vary on the run, but there will typically be 50--150 offpsring (i.e. iterations) per generations. This means that 2500--2.25 &times; 10<sup>4</sup> iterations will be needed to complete the optimization. This can be a **LONG** process!


This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **Help** toolbar button for more details on using R Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
summary(cars)
```

```
##      speed           dist    
##  Min.   : 4.0   Min.   :  2  
##  1st Qu.:12.0   1st Qu.: 26  
##  Median :15.0   Median : 36  
##  Mean   :15.4   Mean   : 43  
##  3rd Qu.:19.0   3rd Qu.: 56  
##  Max.   :25.0   Max.   :120
```


You can also embed plots, for example:


```r
plot(cars)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


