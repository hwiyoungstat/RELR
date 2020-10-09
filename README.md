Robust Extrinsic Local Regression
================

Introduction
------------

This funtion implements the methodology described in the paper

-   *Lee, H.* and Patrangenaru, V. (2020). *Robust Extrinsic Regression Analysis for Manifold Valued Data* [*\[Project Page\]*](https://hwiyoungstat.github.io/RELR.html)

Installation
------------

THis package can be installed with the 'devtools' package:

``` r
library(devtools)
install_github("hwiyoungstat/RELR")
```

Usage
-----

`Ext.med_Shape`, `Multi_RELR` are the primary functions that implement extrinsic median and robust extrinsic local regression on planar shape space, respectively.

``` r
Fit.med <- Ext.med_Shape(X)
```

-   Inputs
    -   `X` : sample data on shape space

``` r
Fit.RELR <- Multi_RELR(grid,X,Response,H)
```

-   Inputs
    -   `grid` : Evaluation points
    -   `X` : Covariate
    -   `Response` : shape (Matrix N\*K)
    -   `H` : bandwidths (column)
