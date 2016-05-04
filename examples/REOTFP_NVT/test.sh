#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Preparing folders
rm -fr output
mkdir output
cd output
mkdir 0 1 2 3 4 5
cd ..

# Run simulation
mpiexec -n 6 namd2 +replicas 6 job0.namd +stdout output/%d/job0.log 2>&1 | tee job0.log
mpiexec -n 6 namd2 +replicas 6 job1.namd +stdout output/%d/job1.log 2>&1 | tee job1.log

# Sorting replicas. For the sort of the replicas, I use sortreplicas code that
# comes with namd distribution. This code expect a name format that has twice
# the number of replica id in the string (in the folder, and in the name), and
# I'm using only 1 (just in the folder). Therefore, this README is a script
# construct hard links to make compatible the names convention.
for j in job1 job0; do
  for i in $(seq 0 5); do
    ln -sf ${j}.dcd output/$i/${j}.${i}.dcd
    ln -sf ${j}.history output/$i/${j}.${i}.history
  done
  ./sortreplicas output/%s/${j} 6 1
  for i in $(seq 0 5); do
    ln -sf ${j}.${i}.sort.dcd output/${i}/${j}_sort.dcd
    ln -sf ${j}.${i}.sort.history output/${i}/${j}_sort.history
  done
done

for file in $(ls output/*/*bsp); do

  f=${file%.*}

  # Decoding the FES
  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA
  
  # Adding x-range
  awk -v N=301 -v min=2.0 -v max=5.0 \
      'END{
        for (i=1;i<=N;i++){
          print min+(max-min)/(N-1)*i,$i
        }
      }' ${f}.LAMBDA > ${f}.fes
       
done

#Computing free energy in equilibrium process of job1
for file in $(ls output/*/job1_ch*.LAMBDA); do

  f=${file%.*}

  awk -v N=301 -v min=2.0 -v max=5.0 \
      '(NR==500){
        for (i=1;i<=N;i++){
          print min+(max-min)/(N-1)*i,$i
        }
      }' $file > ${f}.fes500 

  awk -v N=301 -v min=2.0 -v max=5.0 \
      '(NR==505){
        for (i=1;i<=N;i++){
          print min+(max-min)/(N-1)*i,$i
        }
      }' $file > ${f}.fes505

done
 
for file in $(ls output/*/*dcd); do

  f=${file%.*}

  #Computing traces
  vmd -dispdev text << HERETCL
   
  mol new ../systems/buta.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all
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
gnuplot test0.plt
gnuplot test1.plt

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

