#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Run simulation
namd2 job0.namd > job0.log

# Computing free energy
for file in $(ls *bsp); do

  f=${file%.*}

  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA
  
  sed -n '$s/\(.\)\s\s*\(.\)/\1\n\2/g;$p' ${f}.LAMBDA > ${f}.fes
   
done
 
#Computing the begining of job1
sed -n '1s/\(.\)\s\s*\(.\)/\1\n\2/g;1p' job1_ch0.LAMBDA > job1_ch0.fes1

#Computing traces
for file in $(ls *dcd); do

  f=${file%.*}

vmd -dispdev text << HERETCL
   
  mol new ../butane_files/buta.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all
  mol addfile $file type dcd first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all

  proc distance { i1 i2 } {
    set sel1 [atomselect 0 "index \$i1"]
    set sel2 [atomselect 0 "index \$i2"]
    set coord1 [lindex [\$sel1 get {x y z}] 0]
    set coord2 [lindex [\$sel2 get {x y z}] 0]
    return [vecdist \$coord1 \$coord2]
  } 

  #open file for writing
  set fil [open "${f}.s" w]
  set all [atomselect top all]
  set nframes [molinfo top get numframes]
  puts \$nframes
  for {set frame 0} {\$frame < \$nframes} {incr frame} {
      animate goto \$frame
      puts \$fil [distance 0 10] 
  }
  close \$fil

  exit            

HERETCL

done 

# Compute the plots
gnuplot << 'HEREGNUPLOT'

  min=2.0
  max=4.8
  dr=(max-min)/200

  j=0
  set terminal png;

  set title "CV trace job0"
  set output "job0_".j.".png"; j=j+1
  plot 'job0.s' u ($0+1):1 w l t 'test',\
       'reference/job0.s' u (($0+1)):1 w l t 'reference',\
       '../OTFP/reference/job0.s' u ($0*+1):1 w l t 'without boundaries'  

  set title "Free energy job0"
  set output "job0_".j.".png"; j=j+1
  plot 'job0_ch0.fes' u ($0*dr+min):1 w l t 'test',\
       'reference/job0_ch0.fes' u ($0*dr+min):1 w l t 'reference',\
       '../OTFP/reference/job0_ch0.fes' u ($0*dr+min):1 w l t 'without boundaries' 

HEREGNUPLOT

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

