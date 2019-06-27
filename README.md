# accretionflowproperties
Computes properties of the accretion flow onto a forming gas giant

## Assumptions
In the current version, we assume
1. Spherical symmetry
2. Free-streaming shock temperature everywhere ([Paper II](https://ui.adsabs.harvard.edu/abs/2019arXiv190605869M), Eq. (33))
3. Equal temperatures before and after the shock
4. No dissociation or ionisation

This means that the temperature in the flow could be a bit larger if kappa\*rho\*r is large. Assumption 1 is of course not realistic far from the planet. Assumption 3 has always been seen to hold up to now.

Note that the shock temperature is independent of Assumption 4 but that dissociation or ionisation will change somewhat the temperature structure and the luminosity in the flow.

## Running
To run, call (broken up on two lines just for display purposes!)
```
./accretionflowproperties.awk -v dMdt=<dM/dt (ME/yr)>  -v MP=<Mp (MJ)>  \
    -v RP=<Rp (RJ)>  -v Ld=<Ldownstr (Lsol)> -v Ra=<Racc (RJ)>
```
i.e., for example,
```
./accretionflowproperties.awk -vdMdt=1e-2 -vMP=1 -vRP=2 -vLd=0 -vRa=123
```
The user needs to pick a (reasonable) accretion radius `Ra`. There are a few numerical parameters in the code but nothing very important.

## References
The equations used are to be found in
[Paper I: Marleau et al. (2017)](http://adsabs.harvard.edu/abs/2017ApJ...836..221M) and
[Paper II: Marleau et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019arXiv190605869M)
with the exception of the flow time between two points, which comes from the integral of `dt = dr / v`
with `v = sqrt(2 G MP (1/r - 1/Ra)`. (An appropriate change of variables is needed for the integral...
See e.g. http://www.astro.uu.se/~hoefner/astro/teach/apd_files/apd_collapse.pdf.)

