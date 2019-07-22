_JMJ-V!_

# accretionflowproperties
Computes properties of the accretion flow onto a forming gas giant or its circumplanetary disc (CPD)

_by Gabriel-Dominique Marleau, Universität Tübingen_

## Assumptions
### Planetary-shock case
The file is [accretionflowproperties.awk](accretionflowproperties.awk).

In the current version for the planetary shock case, we assume:
1. Spherical symmetry
2. Free-streaming shock temperature everywhere ([Paper II](https://ui.adsabs.harvard.edu/abs/2019arXiv190605869M), Eq. (33))
3. Equal temperatures before and after the shock
4. No dissociation nor ionisation

This means that the temperature in the flow could be a bit larger if kappa\*rho\*r is large. Assumption 1 is of course not realistic far from the planet. Assumption 3 has always been seen to hold up to now.

Note that the shock temperature is independent of Assumption 4 (except perhaps for rare cases; Marleau et al., in prep.) but that dissociation or ionisation will change somewhat the temperature structure and the luminosity in the flow.

### CPD-shock case
The file is [accretionflowproperties_aboveCPD.awk](accretionflowproperties_above.awk).

In the case of the CPD-surface shock case, because of the very unknown shock geometry, we just set a constant density and temperature above the CPD. Assumptions 2.--3. from above are kept


## Execution
To run, call (shown on two lines just for display purposes)
```
./accretionflowproperties.awk  -v dMdt=<dM/dt (ME/yr)>  -v MP=<Mp (MJ)>  \
    -v RP=<Rp (RJ)>  -v Ld=<Ldownstr (Lsol)>  -v Ra=<Racc (RJ)>
```
i.e., for example,
```
./accretionflowproperties.awk -vdMdt=1e-2 -vMP=1 -vRP=2 -vLd=0 -vRa=123
```
(With `awk`, `-v x=33` or `-vx=33` assigns the value 33 to the variable `x`.) The user needs to pick a (reasonable) accretion radius `Ra`. There are a few parameters in the script. Most are minor but the time zero location makes some difference (try changing it).

Reasonable parameter combinations are maybe
- `Ld = 0--1e-2 # Lsol` (0 is probably usually more relevant)
- `dMdt = 1e-4--1e-2 # ME/an`
- `MP = 1--10 # MJ`
- `RP = 1.5--3 # RJ`


See also [Mordasini et al. (2017), "Characterization of exoplanets from their formation. III. The statistics of planetary luminosities"](http://adsabs.harvard.edu/abs/2017A%26A...608A..72M) for an idea of what internal luminosities could be during formation.

**Reminder**: The output of the script, in particular the temperature, will be approximate (in fact, most likely a lower bound)
because both the formulae (spherical symmetry and even then, simple result only;
in [Paper II](https://ui.adsabs.harvard.edu/abs/2019arXiv190605869M) we give a more exact result)
but also the input (parameter values) are approximate.

## References
The equations used are to be found in
[Paper I: Marleau et al. (2017)](http://adsabs.harvard.edu/abs/2017ApJ...836..221M) and
[Paper II: Marleau et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019arXiv190605869M)
with the exception of the flow time between two points, which comes from the integral of `dt = dr / v`
with `v = sqrt(2 G MP (1/r - 1/Ra))`. (An appropriate change of variables is needed for the integral...
See e.g. http://www.astro.uu.se/~hoefner/astro/teach/apd_files/apd_collapse.pdf.)

