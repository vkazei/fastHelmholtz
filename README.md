# Fast Helmholtz solver
(c) Vladimir Kazei, Oleg Ovcharenko, Dmitry Kabanov (KAUST, 2019)

Solves Helmholtz equation with Enquist, 1977 first (from https://github.com/TristanvanLeeuwen/SimpleFWI) and second order boundary conditions and then compares it with analytic Green's function. 

The discretized equation is solved with a banded matrix solver. Which is the most efficient for shallow, wide models.
![](compare_helm.png)

getA - assemble Helmholtz matrix
getP - project solution to receiver positions
F - froward modeling
testHelmholtz - test on homogeneous model
