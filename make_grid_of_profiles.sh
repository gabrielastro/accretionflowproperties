#!/usr/bin/env bash
# JMJ-V!
# 
# Gabriel-Dominique Marleau, Uni Tübingen
# gabriel.marleau@uni-tuebingen.de
# 
# 
# Driver for a grid of accretion profiles (approximate flow properties)
#   for ./accretionflowproperties.awk
# 
# v.1: 25.10.2019
# 
# ----------------------------------------------------------------------
# 

set -o nounset
set -o pipefail
set -e    # exit on the first mistake

# =======================================================================================================
#  Auxiliary functions
# =======================================================================================================
#
# Radius fit to the population synthesis radii
#   Data from Mordasini et al. (2012, http://adsabs.harvard.edu/abs/2012A%26A...547A.112M)
#   Fits in Aoyama, Marleau, Mordasini, & Ikoma, in prep.
# 
# $1: dM/dt in ME/yr
# $2: MP in MJ
# $3: flag: 1=warm population, 0=cold population
# 
function get_fitted_radius() {
	if [ $# != 3 ]; then
		## do not forget " >&2"!
		echo " *** Wrong number of arguments in get_fitted_radius! ***" >&2
		echo "     Received: $@" >&2
		exit 2
	fi
	
	# calculate the radius with the fit
	awk -v dMdt=$1 -v MP=$2 -v warmPop=$3 'BEGIN{
		lm2 = log(dMdt/1e-2)/log(10);
		if (warmPop==1) {
			## without the backslashes, the rest of the expression is not used...
			printf "%.4f", 0.411-0.244*lm2+3.45*exp(0.762*lm2) \
				+ (-0.489-0.0961*lm2+0.652*exp(0.353*lm2))*(MP-1.0) \
				+ (-0.228-0.00106*lm2+0.226*exp(0.000220*lm2))*(MP-1.0)^2.
		}
		else {
			printf "%.4f", 1.53+0.111*lm2+1.06*exp(0.906*lm2) \
				+ (-0.195-0.0307*lm2+0.0977*exp(0.000695*lm2))*(MP-1.0) \
				+ (-0.250+0.000276*lm2+0.254*exp(0.000214*lm2))*(MP-1.0)^2.
		}
	}'
}

# -------------------------------------------------------------------------------------------------------

# $1: dMdt
# $2: MP
# $3: radius mode or value
function set_radius() {
	if [ $3 == "warm" ]; then
		get_fitted_radius $1 $2 1
	elif [[ ( $3 == "cold" || $3 == "kalt" ) ]]; then
		get_fitted_radius $1 $2 0
	else
		echo $3
	fi
}


# =======================================================================================================
#  Settings (including grid): values and ranges
# =======================================================================================================
# 
# See ./accretionflowproperties.awk for the variables that need to be set
# 
# for each variable: set a value or (in quotes:) a list
#   e.g. M=$(echo 1 {2..10..2})

# the master (summary) file containing the file names and parameter settings:
listofprofilefiles=Ausgabe/Listedatei_Sph+Polar.dat

Ra=200

Ld=0

ffill="1.0 0.3"

# only integer powers for dMdt possible like this; for floating points, change variable to dMdt and adapt below
lgdMdt=$( echo {-5..-1..1} )

#MP=$( echo 1 {2..20..2} )
MP="1 3 5 10 20"

# the radius can be a mode ("cold" or "warm") or a radius, or a set of either or a combination thereof
#   Examples: RP=cold; RP="cold warm"; RP="cold hot 2"
RP="kalt  warm"

# rmax/Ra:
rmaxfact=0.7

# optional settings (currently not loop variables!)
notime=1

take_opacity_into_account=1


# =======================================================================================================
#  Preparations
# =======================================================================================================
# 
# backup the summary file
if [ -f "$listofprofilefiles" ]; then
	mv "$listofprofilefiles" "${listofprofilefiles%.dat}_bak.dat"
fi

# print header; in the loop: only append to the file
#   Make sure the header matches with below!
echo "#-1:profile_file  2:EMPTY 3:EMPTY  4:EMPTY  5:EMPTY  6:EMPTY  \
  7:warm  8:MP  9:lgMPkt  10:Ld     11:ffill   12:dMdt   13:RP \
  14:Ra     15:rmaxfact" > "$listofprofilefiles"

# =======================================================================================================
#  Main loop
# =======================================================================================================
# 
for xRa in $Ra; do
	echo " 1================  Ra = $xRa RJ  ==============="
	for xLd in $Ld; do
		echo "     2================  Ld = $xLd Lsol  ================"
		for xffill in $ffill; do
			echo "       3================  ffill = $xffill  ================"
			for xRP in $RP; do
				echo "          4================  RP = $xRP  ================"
				# set a flag for which radius fit is used
				if [ $xRP == "warm" ]; then warm=1; elif [[ ( $xRP == "cold" || $xRP == "kalt" ) ]]; then warm=0; else warm=manual; fi
				
				for xlgdMdt in $lgdMdt; do
					## currently we loop over the power of dMdt; this is used in the summary file too
					xdMdt=1e$xlgdMdt
					
					echo "         5================  dMdt = $xdMdt ME/yr) ================"
					for xMP in $MP; do
						echo "            6================  MP = $xMP  ================"
					
						# set radius
						theRP=$(set_radius $xdMdt $xMP $xRP)
						echo "            RP="$theRP
						
						# set output file (make sure to creat the subdirectory first)
						outfile=Ausgabe/Profil_$xRP""_fF$xffill"_"$xMP""MJ_MPkt$xdMdt"_"$theRP""RJ_RAkk$xRa""_Ld$xLd.dat
						echo " $outfile"
						
						# calculate the profile
						./accretionflowproperties.awk \
							-v rmaxfact=$rmaxfact -v notime=$notime -v take_opacity_into_account=$take_opacity_into_account \
							-v ffill=$xffill -v dMdt=$xdMdt  -v MP=$xMP  -v RP=$theRP  -v Ld=$xLd -v Ra=$xRa >  $outfile
						
						
						# add this information to the master file
						#  ACHTUNG: xRP is not the same as theRP
						echo -e $outfile '\t' "  x  x  x  x  x  $warm   $xMP   $xlgdMdt   $xLd  $xffill  $xdMdt   $theRP   $xRa   $rmaxfact" \
							>>  "$listofprofilefiles"
					done
					echo
				done
			done
		done
		echo
	done
done

echo " Summary in $listofprofilefiles"
