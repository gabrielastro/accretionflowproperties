#! /usr/bin/awk -f
# JMJ-V!
# 
# Gabriel-Dominique Marleau, Uni Tübingen
# gabriel.marleau@uni-tuebingen.de
# 
# 
# Compute approximate properties of a supersonic accretion flow towards a planet
#   Run this script without arguments to get information
# 
# v.1: 26.06.2019
# v.2: 23.10.2019
# 
# ----------------------------------------------------------------------

# evaporation temperature of the dust from Isella & Natta (2005)
function Tevap(rho) {
	return 2e3*rho^0.0195
}

# 
# mathematical functions
# 
function min(a,b) { return (a<b? a: b) }
function max(a,b) { return (a>b? a: b) }
function log10(x) { return log(x)/log(10.) }

function asin(x) {
    ## thanks to https://stackoverflow.com/questions/1960895/assigning-system-commands-output-to-variable
    K = "perl -E 'use Math::Trig; say asin(" x ")'" 
    K | getline _y   ##  watch out with variable names! Only "x" is a dummy variable
    close(K)
    return _y
}

# the indefinite integral of sqrt(x/(1-x))
##  trick: substitute r/Ra = sin^2(theta) to integrate (e.g.: http://www.astro.uu.se/~hoefner/astro/teach/apd_files/apd_collapse.pdf)
function theIntegral(x) {
    return asin(sqrt(x)) - sqrt(x*(1-x))
}

# set physical and model (up to now: opacity) constants
function setconstants() {
    pi=3.1415926
    G=6.67e-8
    RJ=7.15e9
    MJ=1.898e30
    ME=5.98e27
    an=3.1557e7
    Lsol=3.839e33
    sigSB = 5.67e-5

    # constant opacity for the dust (so uncertain anyway)
    #   per gram gas (i.e., dust-to-gas ratio is built in)
    kappa0 = 3.0  # cm^2/g; cf. Fig. 2 of Paper II

}

function printinfo() {

    print "# ---------------------------------------------------------------------------------------------------"
    print "# Quantities in an approximate 1D accretion profile"
    print "# ---------------------------------------------------------------------------------------------------"
    print "#   Input"
    print "# ---------------------------------------------------------------------------------------------------"
    print "#  ffill: filling factor (-)              = " ffill
    print "#   dMdt: accretion rate (ME/an)          = " dMdt
    print "#     MP: planet mass (MJ)                = " MP
    print "#     RP: planet radius (RJ)              = " RP
    print "#     Ld: downstream luminosity (Lsol)    = " Ld
    print "#     Ra: 'accretion radius' (RJ)         = " Ra
    print "# ---------------------------------------------------------------------------------------------------"
    print "#   Input, more"
    print "#     kappa_nonzero, kappa0  = ", kappa_nonzero, kappa0
    print "#     rmax/Ra, notime, N = ", rmaxfact, notime, N
    print "# ---------------------------------------------------------------------------------------------------"
    print "#   Output"
    print "#     RJ = 7.15e9 cm"
    print "#     t: time to reach the corresponding position with t=t0 (usually =0) at r=r0 (usually =rmax)"
    print "#     t1, t2, t3: estimates of the time to reach the position (just for comparison)"
    print "#     tff: local free-fall time (just for comparison)"
    print "# ---------------------------------------------------------------------------------------------------"
}

BEGIN{
  
  # need constants for printinfo() and main code
  setconstants()
  
  if(length(ffill)*length(dMdt)*length(MP)*length(RP)*length(Ld)*length(Ra) == 0){
    print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "# !!  At least one parameter is not set  !!"
    print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "#"
    print "# ---------------------------------------------------------------------------"
    print "#   Prints out estimated properties of the accretion flow onto a planet      "
    print "# ---------------------------------------------------------------------------"
    print "# "
    print "# References: Marleau et al. (2017, 2019)"
    print "#  http://adsabs.harvard.edu/abs/2017ApJ...836..221M (Paper I)"
    print "#  http://adsabs.harvard.edu/abs/2019ApJ...881..144M (Paper II)"
    print "# "
    print "# Assumptions of the current version:"
    print "#   - Spherical symmetry"
    print "#   - Equal temperatures before and after the shock"
    print "# "
    print "# The model is certainly not realistic at large distances from the planet!"
    print "# "
    print "# Now, in v.2, we do not assume a free-streaming shock temperature everywhere"
    print "#   but rather correct roughly for the effect of the dust opacity."
    print "#   Note that this has not been tested in detail!"
    print "# "
    print "# -------------------"
    print "#   Usage:"
    print "# -------------------"
    print "# "
    print "#  ./accretionflowproperties.awk -v ffill=<ffill> -v dMdt=<dM/dt (ME/yr)>  -v MP=<Mp (MJ)>  -v RP=<Rp (RJ)>  -v Ld=<Ldownstr (Lsol)> -v Ra=<Racc (RJ)>"
    print "# "
    print "# e.g.:   ./accretionflowproperties.awk -vffill=1 -vdMdt=1e-2 -vMP=1 -vRP=2 -vLd=0 -vRa=240 -vrmaxfact=0.7"
    print "#         (-v x=1 or -vx=1  defines and sets a variable)"
    print "# where Ld is the luminosity downstream of the shock = Lint + Lcomp = internal + compression luminosity, and"
    print "#       Ra is the accretion radius (Paper I), which could be computed somewhat self-consistently"
    print "# "
    print "# Optional settings:"
    print "#   -v rmaxfact=<rmax/Ra>    starting point of a parcel of gas as a fraction of the accretion radius [default: 0.99]"
    print "#   -v kappa_nonzero=<0/1>   1 = consider opacity approximately for the temperature profile; 0 = assume free-streaming everywhere [default: 1]"
    print "#   -v notime=<0/1>          1 = skip the computation of the exact time, which is what makes the script so slow [default: 0 = do not skip]"
    print "#   -v N=<number of cells>   number of radial cells between rmax and the planet radius [default: 100]"
    print "# "
    print "# Comments:"
    print "#   * Lint ~ 0 is fine for high accretion rates but is not negligible for lower dM/dt"
    print "#   * Lcomp is not easy to estimate a priori but is also probably not too high..."
    print "#   * With filling factors ffill < 1, one can estimate properties in the accretion column"
    print "#     towards a planet accreting magnetospherically"
    print "#   * Note that for r > (3/4)*Ra, the density reincreases outwards"
    print "# "
    
    printinfo()

    print "# "
    exit 1
  }


# ----------------------------------------------
# 
# optional quantities: further settings
# 

# factor for the outer radius of grid rmax (fraction < Ra)
#   note that the density is monotonic only out to 3/4*Ra;
#   it might be a good idea but is not necessary to start within this radius.
if (length(rmaxfact) == 0) { rmaxfact = 0.99 }

# take opacity into account? By default, yes (this is more accurate)
#   This can increase the temperature in the flow
if (length(kappa_nonzero) == 0) { kappa_nonzero = 1 }

# skip computation of exact time? This computation makes the script slow
#   default: do compute
if (length(notime) == 0) { notime = 0 }

# number of radial points minus 1
if (length(N) == 0) { N = 100 }

# ----------------------------------------------

# print information before unit conversion of variables
printinfo()

# ----------------
#  convert units of input and of some derived variables
# ----------------
dMdt  = dMdt * ME/an
MP    = MP   * MJ
RP    = RP   * RJ
Ld    = Ld   * Lsol
Ra    = Ra   * RJ
ffill = ffill  # this is obviously dimensionless

# outer radius
rmax = rmaxfact * Ra

# luminosity just above the shock
Ltot = G*MP*dMdt/RP*sqrt(1-RP/1./Ra) + Ld

# free-streaming temperature at the shock (downstream + accretion)
Tshff = ( Ltot/(4*pi*RP^2.*ffill) / (4*sigSB) )^0.25

# check grid
if(N / log10(Ra/RP) < 10){
    print "# * ACTHUNG: " N " points for " log10(Ra/RP) " dex: not very good resolution for integration! *"
}

# print more information
print "#   More information"
print "#     Ltot (Lsol), Tshff (K): ", Ltot/Lsol, Tshff

# 
# print out every layer from outside in, from rmax down to RP
#   apart for the numerical integration of the time,
#   every layer is independent of the others
#  
for(i=0;i<=N;i++){
  
  # radial position
  #   RP <= r <= rmax < Ra
  ## save for integration
  r_old = r
  ## 
  ## logarithmic grid
  r = rmax *(RP/rmax)^( i/(N*1.) )
  ## 
  ## linear grid
  #r = rmax - i*(rmax-RP)/(N*1.)
  
  # time -- t is exact and the others are approximations
  ## 
  ## define time zero to be when a parcel of gas is at radius rmax
  ##   there, it has a velocity v0 = vff(r0, Ra) = sqrt(2 G MP (1/r - 1/Ra))
  t0 = 0
  r0 = rmax
  ## 
  ## variables which will show up many times below
  x = sqrt(Ra^3/(2*G*MP))
  y  = r/Ra
  y0 = r0/Ra
  ##
  ## free-fall from Ra, exact, analytic
  ##  calling theIntegral() is what makes the script slow
  if(notime==0) {
    t  = t0 + x*( theIntegral(y0) - theIntegral(y) )
  } else {
    t  = -33.  # serves as a flag
  }
  ## 
  ## free-fall from infinity vff = -sqrt(2 G MP / r), exact result, analytic
  ##   not consistent with the rest (Ra = infty instead of input value)
  t1 = t0 + x* 2/3.*( y0^1.5 - y^1.5 )
  ## 
  ## free-fall from Ra, approximate result (Taylor expansion, for r << Ra), analytic
  ##   consistent vff = -sqrt(2 G MP (1/r - 1/r0) )
  t2 = t0 + x*( 2/3.*( y0^1.5 - y^1.5 ) + 1/5.*( y0^2.5 - y^2.5 ) )
  ## 
  ## free-fall from Ra, exact formula, simple numerical integration
  if(i==0){
    # initialise
    t3 = t0
  } else {
    # dt = dr / -v
    #   use logarithmic middle r to evaluate v
    rmid = sqrt(r*r_old)
    t3 = t3 + (r_old - r) / sqrt(2*G*MP*(1./rmid - 1./Ra))
  }
  ##
  ## for comparison only, the "free-fall time" at each radius
  ##   the mass of the gas in the flow is negligible
  tff = sqrt(3*pi/(32*G* MP/(4*pi/3.*r^3.)))
  
  # free-fall (=no pressure gradient) velocity,
  #   not necessarily from infinity
  vff = - sqrt(2*G*MP*(1./r - 1./Ra))
  
  # free-fall density
  rhoff = dMdt/(4*pi*r^2.*ffill * -vff)
  
  # temperature, assuming T \propto (r/RP)^{-1/2} and a free-streaming shock temperature
  T = Tshff * (RP/r)^0.5
  
  # 
  # multiply by opacity factor (Eq. 32 of Paper II)
  ##   (compare with Fig. 9b of Paper II)
  ##   Note: instead of finding the temperature implicitly self-consistently,
  ##   we use a rough approach
  if(kappa_nonzero==1){
    
    Tdestloc = Tevap(rhoff)
    # use a step function for the opacity with a height ~ typical value
    kappa = (T>=Tdestloc? 0 : kappa0)
    xkap = kappa*rhoff*r
    
    # from Eq. (35) of Paper II:
    ##   for the Ensman (1994) flux limiter:
    #fact = (1 + 3./2.*xkap)^0.25
    ##   for the Levermore & Pomraning (1981) flux limiter:
    fact = ( (1 + 3./2.*xkap + 3./2.*xkap^2.)/(1+xkap) )^0.25
    
    # this is the rough fit: below the destruction temperature,
    ##   take the minimum (with smoothing) of the two
    if(T <= Tdestloc) {
      #T = min (Tdestloc, T * fact)
      ## a strong power matches well at least for the case of Fig. 9 in Paper II
      p = 10.
      T = ( 1./Tdestloc^p + 1/(T*fact)^p )^(-1/p)
    }
  
  }  # end of kappa_nonzero
  
  if(i==0) {
    # print header line
    print "# ---------------------------------------------------------------------------------------------------"
    print "#-1:i      2:t/s    3:r/RJ   4:rho/(g/cm^3)     5:T(K)  6:v/(km/s)       7:t1/s      8:t2/s      9:t3/s    10:tff/s"
  }
  
  # print the layer
  printf " %4d  %10.3e  %8.3f  %15.4e  %9.2f  %10.3f   %10.3e  %10.3e  %10.3e  %10.3e\n",
      i, t, r/RJ, rhoff, T, vff/1e5, t1, t2, t3, tff
  
}  # end of main loop

}
# done!
