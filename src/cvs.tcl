
# load the C-module
load ${CFACV_BASEDIR}/cvs.so cvs

proc parse_pdb {pdbfile pdbarray} {
  upvar $pdbarray pdb 

  # Read pdb and use it to construct pdbarray

  # Open and read the pdb. Put all it in lines variable.
  set inStream [open $pdbfile r]
  set data [read $inStream]
  close $inStream
  set lines [split $data \n]

  set j 0 

  foreach line $lines {

    # Check if we have an ATOM. NAMD does not pay atention to pdb atom index so
    # we will just count.
    set cardstr [string range $line 0 3]
    if {$cardstr ne "ATOM"} continue
 
    # Retrive atom index
    incr j
 
    # Check and retrive center index
    set i [expr int([string trim [string range $line 60 65]])]
    if {$i == 0} continue

    # Save center data
    if {![info exists pdb($i.atoms)]} {incr pdb(nctrs)}
    lappend pdb($i.atoms) $j

    # Save atom data
    incr pdb(natoms)
    set pdb(a$j.resname) [string trim [string range $line 17 19]]
    set pdb(a$j.x)       [string trim [string range $line 30 37]]
    set pdb(a$j.y)       [string trim [string range $line 38 45]]
    set pdb(a$j.z)       [string trim [string range $line 46 53]]
    set pdb(a$j.mass)    [string trim [string range $line 54 59]] ; #occupancy
    # set pdb(occustr) [string trim [string range $line 54 59]]

    # Save the pdb line
    lappend pdbline $line
  }
}
  
proc cvs_setup { pdbarray } {
  global ds
  global cv
  upvar $pdbarray pdb

  for {set i 1} {$i<=$cv(num)} {incr i} {


    # Defaults and consistency in boundaries
    if {![info exists cv($i.bound)]} {set cv($i.bound) NADA}
    set cv($i.bound) [string toupper $cv($i.bound)]
    switch $cv($i.bound) {
      NADA {}
      SOFTUPPER {
        if {![info exists cv($i.max)]}    {error "CFACV) ERROR: cv need min value"}
        if {![info exists cv($i.boundk)]} {error "CFACV) ERROR: cv need boundk value"}
      }     
      SOFTLOWER {
        if {![info exists cv($i.min)]}    {error "CFACV) ERROR: cv need max value"}
        if {![info exists cv($i.boundk)]} {error "CFACV) ERROR: cv need boundk value"}
      }     
      SOFT {
        if {![info exists cv($i.min)]}    {error "CFACV) ERROR: cv need max value"}
        if {![info exists cv($i.max)]}    {error "CFACV) ERROR: cv need min value"}
        if {![info exists cv($i.boundk)]} {error "CFACV) ERROR: cv need boundk value"}
      }
      default {error "ERROR: $cv($i.bound) is not a valid boundary type."} 
    }

    if {![info exists cv($i.boundk)]} {set cv($i.boundk) 0.}
    if {![info exists cv($i.min)]}    {set cv($i.min) 0.}
    if {![info exists cv($i.max)]}    {set cv($i.max) 0.}
          
     
    # Fill centers index 0-based list needed for C
    set ind {}
    foreach ctr $cv($i.centers) {
      if {$ctr > $pdb(nctrs)} {error "index number higher than the number of centers ($pdb(ncntr))"}
      if {$ctr < 1} {error "index number should be grater than 0"}
      lappend ind [expr $ctr-1]
    }
    set nInd [llength $ind]
    set pind [intListToArray $ind]
    

    # Output
    if {[info exists cv($i.outfile)]} {
      if {![info exists cv($i.outfreq)]}  {error "CFACV) ERROR: cv need output frequency"}
    } else {
      if {[info exists cv($i.outfile)]}  {error "CFACV) ERROR: cv need output file"}
      set cv($i.outfile) ""
      set cv($i.outfreq) -1
    }
    

    # The CV type
    set typ [lindex $cv($i.type) 0] 


    # Allocate the cv
    set cv($i.address) [DataSpace_add_cv $ds $typ $nInd $pind $cv($i.min) $cv($i.max) $cv($i.bound) $cv($i.boundk) $cv($i.outfile) $cv($i.outfreq)]


    # Post process of the cv structure
    switch $typ {
      ZSDCIRCLE {
        set_zsd_circle [lindex $cv($i.type) 1] [lindex $cv($i.type) 2] [lindex $cv($i.type) 3] [lindex $cv($i.type) 4]
      }
      ZSDXRANGE {
        set_zsd_circle [lindex $cv($i.type) 1] 0.0 [lindex $cv($i.type) 2] [lindex $cv($i.type) 3]
      }
      ZSDRING {
        set_zsd_ring [lindex $cv($i.type) 1] [lindex $cv($i.type) 2] [lindex $cv($i.type) 3]  [lindex $cv($i.type) 4] [lindex $cv($i.type) 5]
      }
      RMSD {
        cv_readref $i [lindex $cv($i.type) 1]
      }
      LINE {
        cv_readref $i [lindex $cv($i.type) 1]
        cv_readref2 $i [lindex $cv($i.type) 2]
       set_line $cv($i.address)
      }
    }
    
  }
}


proc cv_readref { i pdbfile} {
  #Read the pdb reference and set the corresponding array in the cv structure.
  global cv
   
  parse_pdb $pdbfile pdb

  set k 0
  foreach j $cv($i.centers) {

    # compute the center of mass of the center
    set cmpos [get_cmpos $pdb($j.atoms) pdb]

	  ListToArray_Data [cv_access_ref $cv($i.address) $k] $cmpos
    
    incr k
  }         
    
}
 
proc cv_readref2 { i pdbfile } {
  #Read the pdb reference and set the corresponding array in the cv structure.
  global cv
   
  parse_pdb $pdbfile pdb

  set k 0
  foreach j $cv($i.centers) {

    # compute the center of mass of the center
    set cmpos [get_cmpos $pdb($j.atoms) pdb]

	  ListToArray_Data [cv_access_ref2 $cv($i.address) $k] $cmpos
    
    incr k
  }         
    
}
  
proc get_cmpos { atomlist pdbarray } {
  # compute the center of mass of the atomlist
  # TODO: make this a c routine. Might be construct a structure with pdbarray.
  upvar $pdbarray pdb 
 
  set mass 0.
  set xpos 0.
  set ypos 0.
  set zpos 0.

  foreach i $atomlist {

    set mass [expr $mass+$pdb(a$i.mass)]
    set xpos [expr $xpos+$pdb(a$i.x)*$pdb(a$i.mass)]
    set ypos [expr $ypos+$pdb(a$i.y)*$pdb(a$i.mass)]
    set zpos [expr $zpos+$pdb(a$i.z)*$pdb(a$i.mass)]
                          
  }

  return [list [expr $xpos/$mass] [expr $ypos/$mass] [expr $zpos/$mass]]
}
