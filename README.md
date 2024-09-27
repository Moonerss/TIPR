
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TIPR

<!-- badges: start -->
<!-- badges: end -->

The goal of TIPR is to run the analysis of TIP web server.

## Installation

You can install the development version of TIPR like so:

``` r
remotes::install_github('Moonerss/TIPR')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TIPR)

data('expression')

input_expr <- log2(expression + 1)
res <- TIP(expression = input_expr)
#> ℹ Filtering gene set ...
#> ℹ Random 100 times ...
#> ℹ Calculate activity score ...
#> ℹ Normalize activity score ...
#> ✔ Done
head(res)
#>                                  Sample1      Sample2      Sample3      Sample4
#> Step1                           9.950372     9.950372     9.950372     9.950372
#> Step2                        -421.618478 -1081.840444  -755.772366 -1310.635579
#> Step3                       -2869.732823  -786.838114  -822.262994  -469.117654
#> Step4.B cell.recruiting      5229.557410   114.195986  3958.385906  2617.812689
#> Step4.Basophil.recruiting   -3825.428263 -4217.397148 -2429.917533 -1576.612513
#> Step4.CD4 T cell.recruiting  7732.909107  7346.025678  5091.805751  4134.347672
#>                                  Sample5
#> Step1                           9.950372
#> Step2                       -1706.976618
#> Step3                        -475.540093
#> Step4.B cell.recruiting      1145.479180
#> Step4.Basophil.recruiting    -604.034790
#> Step4.CD4 T cell.recruiting  3978.223482
```
