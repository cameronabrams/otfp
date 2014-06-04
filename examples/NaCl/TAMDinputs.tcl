#!/usr/local/bin/vmd -e

set istru NaCl15_wb30

mol load psf ${istru}.psf pdb ${istru}.pdb     


# empty board
set all [atomselect top all]     
$all set beta 0          
$all set occupancy 0

# selection
set sel [atomselect top "ions"]

# Mass
$sel set occupancy [$sel get mass]    

# Beta
set n [$sel num]
set l {}
for {set i 1} {$i <= $n} {incr i} { lappend l $i.00 }
$sel set beta $l

# pdb
$all writepdb label.pdb

# cv
set ofp [open "cv.inp" "w"]
for {set i 0} {$i < $n} {incr i} {
  puts $ofp "CARTESIAN_X $i"
  puts $ofp "CARTESIAN_Y $i"
  puts $ofp "CARTESIAN_Z $i"
}




exit     

