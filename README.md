## Generalized Partially Linear Single Index Model with  Measurement Error, Instruments and Binary Response
### Guangren Yang, Qianqian Wang, Xia Cui and Yanyuan Ma

* **Abstract:** 
Partially linear generalized single index models are widely used and have attracted much attention in the literature. However, when the covariates are subject to measurement error, the problem is much less studied. On the other hand, instrumental variables are important elements in studying many errors-in-variables problems. We use the relation between the unobservable variables and the instruments to devise consistent estimators for partially linear generalized single index models with binary response. We establish the consistency, asymptotic normality of the estimator and illustrate the numerical performance of the method through simulation studies and a data example. 


* **How to use**
File ending in .f90 contains code for the simulations.
R file contains code for plots and calculated results.
The gfortransubs file contains functions used in Fortran.

  To compile please use: 

  ```
  gfortran -O2 Fortran_code.f90 ~/gfortransubs/* -o simu1
  ```
