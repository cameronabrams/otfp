#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

getnamdset () {
  sed -n 's/^ *set  *'$1' .*/\0\nputs $'$1'/p' $2 | tclsh
}
 
# # Run simulation
# namd2 job0.namd 2>&1 | tee job0.log

# Run simulation
# namd2 job1.namd 2>&1 | tee job1.log

# Computing free energy.
for file in $(ls *bsp); do

  f=${file%.*}

  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA

  awk  -v line=$(wc -l < ${f}.LAMBDA) -v j=0 \
     -v Nx=151 -v Ny=151 \
     -v dx=0.0418879020478639 -v dy=0.0418879020478639 \
     -v xmin=-3.14159265358979 -v ymin=-3.14159265358979 \
     '(line==NR){
        for (j=0;j<Ny;j++){
          for (i=0;i<Nx;i++){
            a=i+j*Nx+1
            print ymin+dy*j,xmin+dx*i,$a
          }
          print ""
        }
      }' ${f}.LAMBDA > ${f}.fes
done

# #Computing the begining of job1
# sed -n '1s/\(.\)\s\s*\(.\)/\1\n\2/g;1p' job1_ch0.LAMBDA \
#   | awk '{print (NR-1)*'$dr'+'$rmin',$1}' \
#   > job1_ch0.fes1

#Computing traces
for file in $(ls *dcd); do

  f=${file%.*}

vmd -dispdev text << HERETCL
   
  mol new ../systems/ala2.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all
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
      puts \$fil "[measure dihed {2 0 16 14}] [measure dihed {0 16 14 13}]"
  }
  close \$fil

  exit            

HERETCL

done 

# Compute the plots
gnuplot test.plt

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

