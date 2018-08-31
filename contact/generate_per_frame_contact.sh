#!/bin/bash

if test $# -lt 1; then

  echo "Usage: $0 mdfiles_list [ncore]"
  exit 1

fi
ulimit -s unlimited

################################################
##    USER CAN CHANGE FOLLOWING PARAMETERS    ##
##         TO CUSTOMIZE THIS SCRIPT           ##
##                                            ##
## NOTE: DO NOT PUT SPACE BEFORE OR AFTER "=" ##
################################################ 

# Distance cutoff to define contact.
cutoff=4.5    

# The program to do the calculation.
prog=./generate_per_frame_contact_v4_ncdump

# Atoms you want to strip in addition to hydrogen,
# such as substrate.
# Follow Amber atom mask syntax.
strip2="" 
#strip2=":1-11"

# Clean up temporary files after calculation?
cleanup=yes   

# Keep stipped trajectory files?
keeptraj=no


###################################
## DO NOT CHANGE FOLLOWING LINES ##
###################################

top=$( sed -n '1'p $1 )
files=( $( sed -n '2,$'p $1 ) )
if test $# -ge 2; then
  ncore=$2
else
  ncore=1
fi

if test ! -x $prog; then
  echo Program not found.
  exit 1
fi

if test -z $( which cpptraj 2>/dev/null ); then
  echo cpptraj not found.
  exit 1
fi 

## - generate topology pdb
echo "trajin $files 1 1" > __tmp_cpptraj0.in
echo "strip @H*" >> __tmp_cpptraj0.in
if test ! -z $strip2; then 
  echo "strip $strip2" >> __tmp_cpptraj0.in
fi
echo "center !:WAT,Na+,Cl-,K+" >> __tmp_cpptraj0.in
echo "image center familiar" >> __tmp_cpptraj0.in
echo "trajout __tmp_top.pdb pdb nobox" >> __tmp_cpptraj0.in
cpptraj $top < __tmp_cpptraj0.in > /dev/null

## - renumber residues in the pdb file
grep '^ATOM' __tmp_top.pdb | cut -c-22 > __tmp_top0.pdb
grep '^ATOM' __tmp_top.pdb | cut -c23-26 | awk 'NR==1{a=$1} {printf "%4d\n", $1-a+1}' > __tmp_top1.pdb
grep '^ATOM' __tmp_top.pdb | cut -c27- | paste -d '' __tmp_top0.pdb __tmp_top1.pdb - > __tmp_top2.pdb
mv __tmp_top2.pdb __tmp_top.pdb


##- if stripped trajecotry files exist, check if we 
##- can use them or not.
bGenTraj=yes
trajs=trajs_noh
if test -d $trajs; then
  if test $( ls $trajs/* | wc -l ) -ne $ncore; then
    echo Directory $trajs exists but does not match ncore.
    read -p "Do you want to re-do trajectory stripping? (yes/no)" ans
    if test $ans == "yes"; then
      rm -rf $trajs
    elif test $ans == "no"; then
      bGenTraj=no
      ncore=$( ls $trajs/* | wc -l )
    else
      exit 1
    fi
  else
    echo Directory $trajs exists and will be used for calculation.
    bGenTraj=no
  fi
fi 
  
if [ $bGenTraj == "yes" ] && [ -r $trajs ]; then
  echo File $trajs exists. Please rename/remove it to continue.
  exit 1
fi

##- how many frames?
echo ${files[@]} | xargs -n 1 ncdump -h > __tmp_heads
nframes=$( grep 'frame\s*=' __tmp_heads | sed 's/.*(\s*\([0-9]*\).*/\1/' |\
   awk 'BEGIN{a=0}{a=a+$1}END{print a}' )

##- how many frames for each core? 
load=( $( yes $(( $nframes / $ncore )) | head -n $ncore ) )
mod=$(( $nframes % $ncore ))
for((i=0; i<$mod; i++)); do
  let load[i]=${load[i]}+1
done



## 1. Strip hydrogen atoms from trajectories.
if test $bGenTraj == "yes"; then

echo "Processing trajectory files ($nframes frames) ..."

mkdir $trajs
rm -f __tmp_cpptraj.in
for i in ${files[@]}; do
  echo "trajin $i" >> __tmp_cpptraj.in
done
echo "strip @H*" >> __tmp_cpptraj.in
if test ! -z $strip2; then 
  echo "strip $strip2" >> __tmp_cpptraj.in
fi
echo "center !:WAT,Na+,Cl-,K+" >> __tmp_cpptraj.in
echo "image center familiar" >> __tmp_cpptraj.in
#echo "trajout __tmp_top.pdb pdb nobox start 1 stop 1" >> __tmp_cpptraj.in

j=0
for((i=0; i<${#load[@]}; i++)); do
  jj=$( echo $i | awk '{printf "%03d", $1+1}' )
  echo -n "trajout $trajs/traj${jj}.nc netcdf " >> __tmp_cpptraj.in
  echo -n "start $((j+1)) stop $(($j+${load[i]})) " >> __tmp_cpptraj.in 
  echo    "nobox noforce novelocity notemperature notime" >> __tmp_cpptraj.in
  let j=j+${load[i]}
done

s0=$( date +"%s" )

cpptraj $top < __tmp_cpptraj.in > /dev/null

if test $? -ne 0; then
  echo Error in running cpptraj - Check paths to MD files.
  exit 1
fi

s1=$( date +"%s" )
echo Finished in $(( $s1 - $s0 )) seconds.
echo

fi

## 2. Generate contact trajectory
echo Generating contact trajectory... 

##- get basic information about the system.
nres=$( grep '^ATOM' __tmp_top.pdb | cut -c 23-26 | uniq | wc -l )
natom=$( grep '^ATOM' __tmp_top.pdb | wc -l )

tsf=$( date +"%s" )
dynamics=__tmp_dynamics_$tsf
mkdir $dynamics

s0=$( date +"%s" )
for((i=0; i<$ncore; i++)); do
  jj=`echo $i | awk '{printf "%03d", $1+1}'`
  ofile="$dynamics/dynamics$jj"
  infile="$trajs/traj${jj}.nc"
  ncdump -v coordinates $infile |  sed -e '1,/coordinates =/'d -e '/}/'d -e 's/[,;]//g' |\
    $prog __tmp_top.pdb $nres $natom ${load[i]} $cutoff $ofile &
done

trap "kill -15 -$$; exit 1" TERM INT KILL

wait

if test $? -ne 0; then
  echo Error in running $prog.
  exit 1
fi

s1=$( date +"%s" )
echo Finished in $(( $s1 - $s0 )) seconds.
echo

##- collect results.
cat $dynamics/* > contact.trajectory


## 3. clean up
if test $cleanup == "yes"; then
  echo Cleaning up...
  rm -rf __tmp_heads __tmp_cpptraj0.in __tmp_cpptraj.in \
         __tmp_top*.pdb $dynamics
  if test $keeptraj != "yes"; then
    rm -rf $trajs
  fi
fi

echo done.
