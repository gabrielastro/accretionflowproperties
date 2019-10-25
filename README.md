_JMJ-V!_

# accretionflowproperties
Computes approximate properties of the accretion flow onto a forming gas giant:
density, temperature, velocity as a function of radius and time since the beginning
of the fall.

_by Gabriel-Dominique Marleau, Universität Tübingen_

The main file is [accretionflowproperties.awk](accretionflowproperties.awk).

## Assumptions

In the current version for the planetary shock case, we assume:
1. Spherical symmetry
2. Equal temperatures before and after the shock
3. No dissociation nor ionisation

- Assumption 1 is of course not realistic far from the planet. In the case of ffill<1 (see below), the profile should be thought of as the average profile over the accreting solid angle.
- Assumption 2 has always been seen to hold up to now.
- The shock temperature is independent of Assumption 3 (except perhaps for rare cases; Marleau et al., in prep.) but dissociation or ionisation will change somewhat the temperature structure and the luminosity in the flow.

In the current version (v.2), we now take roughly into account the effect of the dust opacity, which is important if locally kappa\*rho\*r is >~ 1; this will typically be the case below the dust destruction tempeature (see [Paper II](http://adsabs.harvard.edu/abs/2019ApJ...881..144M), e.g., Fig. 7). Given the uncertainties in the dust model, we use a constant kappa = 3 cm^2/g (cross-section per gram of _gas_). Caveat emptor: the simple smoothing used here (instead of actual solving implicity) has not been tested extensively!


## Input and output
### Required input
The required input quantities are:
- filling factor (fraction of planet surface covered by accretion) ffill
- accretion rate dM/dt
- planet mass MP
- planet radius = shock position RP
- downstream luminosity Ld
- accretion radius Ra

Ld is the luminosity left of the shock = internal luminosity + post-shock compression of the material,
and needs to be determined by some other means; taking Ld=0 (in the sense << L_acc) is usually fine.
See also [Berardo et al. (2017)](http://adsabs.harvard.edu/abs/2017ApJ...834..149B).

Ra is the "accretion radius" in the sense that v_free-fall(r) = sqrt[ 2 G MP * (1/r - 1/Ra) ].


### Optional settings
- rmaxfact=<rmax/Ra>    starting point of a parcel of gas as a fraction of the accretion radius; default: 0.99
- kappa_nonzero=<0/1>   1 = consider opacity approximately for the temperature profile; 0 = assume free-streaming everywhere; default: 1
- notime=<0/1>          1 = skip the computation of the exact time, which is what makes the script so slow; default: 0 = do not skip
- N=<number of cells>   number of radial cells between rmax and the planet radius; default: 100

Since the longest part of the infall is the early phase, the "time zero" location (`rmaxfact`)
makes some difference on the total time to reach the planet (try changing it).
Note that all other quantities are independent of time, however.


### Output
The output quantities are, as a function of time or position:
- time since the accretion radius (this makes the script slow; can be deactivated)
- radial position (distance from planet centre)
- density
- temperature
- velocity
- crude estimates of the time coordinate, only for comparison
- local free-fall time

The output will be printed to screen.

## Execution
To run, call (shown on two lines just for display purposes)
```
./accretionflowproperties.awk  -v ffill=<ffill>  -v dMdt=<dM/dt (ME/yr)>  -v MP=<Mp (MJ)>  \
    -v RP=<Rp (RJ)>  -v Ld=<Ldownstr (Lsol)>  -v Ra=<Racc (RJ)>
```
i.e., for example,
```
./accretionflowproperties.awk -vffill=1 -vdMdt=1e-2 -vMP=1 -vRP=2 -vLd=0 -vRa=240 -v rmaxfact=0.7
```
(With `awk`, `-v x=33` or `-vx=33` assigns the value 33 to the variable `x`.)

- Note:
The user needs to pick a (reasonable) accretion radius `Ra`, e.g. Ra = 1/3 of the Hill sphere RH;
this is 240 RJ for Jupiter at 5 au and 650 RJ for 20 MJ at 5 au.

Reasonable parameter ranges are maybe:
- `ffill = 0.01--1` (1: spherically symmetric; e.g., 0.01: magnetospheric accretion)
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
[Paper II: Marleau et al. (2019)](http://adsabs.harvard.edu/abs/2019ApJ...881..144M)
with the exception of the flow time between two points, which comes from the integral of `dt = dr / v`
with `v = sqrt(2 G MP (1/r - 1/Ra))`. (An appropriate change of variables is needed for the integral...
See e.g. http://www.astro.uu.se/~hoefner/astro/teach/apd_files/apd_collapse.pdf.)


## Contact
Corrections, comments, requests, and suggestions are all welcome! Please contact me at `uni-tuebingen.de`, with `gabriel.marleau` in front.
