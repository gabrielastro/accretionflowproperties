#! /usr/bin/awk -f
# JMJ-V!
# 
# 26.06.2019 Gabriel-Dominique Marleau, Uni Tübingen
# gabriel.marleau@uni-tuebingen.de
# 
# Run without arguments to get information
# 

function log10(x) {
    return log(x)/log(10.)
}

# arcsine
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

function header() {

    print " # Simple accretion profile, free-streaming"
    print " # -------------------"
    print " #   Input"
    print " # -------------------"
    print " #  ffill: filling factor (-)              = " ffill
    print " #   dMdt: accretion rate (ME/an)          = " dMdt
    print " #     MP: planet mass (MJ)                = " MP
    print " #     RP: planet radius (RJ)              = " RP
    print " #     Ld: downstream luminosity (Lsol)    = " Ld
    print " #     Ra: 'accretion radius' (RJ)         = " Ra
    print " # Numerical: notime, N, rmax/7.15e9: ", notime, N, rmax/7.15e9
    print " # -------------------"
    print " #   Output"
    print " #     RJ = 7.15e9 cm"
    print " #     t: time to reach the corresponding position with t=t0 (usually =0) at r=r0 (usually =rmax)"
    print " #     t1, t2, t3: estimates of the time to reach the position (just for comparison)"
    print " #     tff: local free-fall time (just for comparison)"
    print " # -------------------"
    print " #-1:i      2:t/s    3:r/RJ   4:rho/(g/cm^3)     5:T(K)  6:v/(km/s)       7:t1/s      8:t2/s      9:t3/s    10:tff/s"
}

BEGIN{
  if(length(dMdt)*length(MP)*length(RP)*length(Ld)*length(Ra)*length(ffill) == 0){
    print " # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print " # !!  At least one parameter is not set  !!"
    print " # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print " #"
    print " # ---------------------------------------------------------------------------"
    print " #   Prints out estimated properties of the accretion flow onto a planet      "
    print " # ---------------------------------------------------------------------------"
    print " # "
    print " # References: Marleau et al. (2017, 2019)"
    print " #  http://adsabs.harvard.edu/abs/2017ApJ...836..221M (Paper I)"
    print " #  https://ui.adsabs.harvard.edu/abs/2019arXiv190605869M (Paper II)"
    print " # "
    print " # Assumptions of the current version:"
    print " #   - Spherical symmetry"
    print " #   - Free-streaming shock temperature everywhere (Paper II, Eq. (33))"
    print " #   - Equal temperatures before and after the shock"
    print " # "
    print " # The model is certainly not realistic at large distances from the planet!"
    print " # "
    print " # -------------------"
    print " #   Usage:"
    print " # -------------------"
    print " # "
    print " #  ./accretionflowproperties.awk -v ffill=<ffill> -v dMdt=<dM/dt (ME/yr)>  -v MP=<Mp (MJ)>  -v RP=<Rp (RJ)> \ "
    print " #           -v Ld=<Ldownstr (Lsol)> -v Ra=<Racc (RJ)>"
    print " # "
    print " # e.g.:   ./accretionflowproperties.awk -vffill=1 -vdMdt=1e-2 -vMP=1 -vRP=2 -vLd=0 -vRa=123"
    print " #         (-v x=1 or -vx=1  defines and sets a variable)"
    print " # where Ld is the luminosity downstream of the shock = Lint + Lcomp = internal + compression luminosity, and"
    print " #       Ra is the accretion radius (Paper I), which could be computed somewhat self-consistently"
    print " # "
    print " # Use the option `-v notime=1` to skip the computation of the exact time, which is what makes the script so slow"
    print " # "
    print " # Comments:"
    print " #   * Lint ~ 0 is fine for high accretion rates but is not negligible for lower dM/dt"
    print " #   * Lcomp is not easy to estimate a priori but is also probably not too high..."
    print " #   * With filling factors ffill < 1, one can estimate properties in the accretion column"
    print " #     towards a planet accreting magnetospherically"
    print " #   * Note that for r > (3/4)*Ra, the density reincreases outwards"
    print " # "

    header()

    print " # "
    exit 1
  }


# 
# optional quantities (numerical settings)
# 
# number of radial points minus 1
if (length(N) == 0) { N = 100 }

# skip computation of exact time? This computation makes the script slow
#   default: do compute
if (length(notime) == 0) { notime = 0 }


# print information before unit conversion
header()

# ----------------
#  constants
# ----------------
pi=3.1415926
G=6.67e-8
RJ=7.15e9
MJ=1.898e30
ME=5.98e27
an=3.1557e7
Lsol=3.839e33
sigSB = 5.67e-5

# ----------------
#  convert units
# ----------------
dMdt = dMdt*ME/an
MP = MP * MJ
RP = RP * RJ
Ld = Ld * Lsol
Ra = Ra * RJ
ffill = ffill  # this is obviously dimensionless


# outer radius of grid, < Ra
#   note that the density is monotonic only out to 3/4*Racc but this does not matter
rmax = 0.99 *Ra

# luminosity just above the shock
Ltot = G*MP*dMdt/RP + Ld

# free-streaming temperature at the shock (downstream + accretion)
Tshff = ( Ltot/(4*pi*RP^2.*ffill) / (4*sigSB) )^0.25

# check grid
if(N / log10(Ra/RP) < 10){
    print " # * ACTHUNG: " N " points for " log10(Ra/RP) " dex: not very good resolution for integration! *"
}

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
  ##  calling theIntegral() is what makes the script slow!
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
  
  # print to screen
  printf " %4d  %10.3e  %8.3f  %15.4e  %9.2f  %10.3f   %10.3e  %10.3e  %10.3e  %10.3e\n",
      i, t, r/RJ, rhoff, T, vff/1e5, t1, t2, t3, tff
  
}  # end of main loop

}
# done!
