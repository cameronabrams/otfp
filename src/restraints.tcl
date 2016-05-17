# A "restraint" is the value of a linear combination of collective
# variables:
# 
# R_l = \sum_{k=1}^{nCV} C_{kl} \Theta_k
#
# where \Theta_k is the value of the k'th collective variable, and 
# C_{kl} is the coefficient multipliying \Theta_k in the l'th
# restraint.
#
# A restraint measure is 
# 
# Z_l = R_l - B_1
# 
# where B_1 is the target value of the restraint.
#
# The potential energy of the restraint is harmonic:
#
# U_l = (1/2) k_l Z_l^2
#
# so the force on a degree of freedom x is
#
# f = -\grad_x U_l = -k_l Z_l \grad_x Z_l
#
#
# Example:  If collective variable \Theta_0 is a bond between centers 0 and 1
# that we wish to restrain to a distance of B_0 = 10 A using restraint R_0, and there
# are a total of 3 collective variables defined in the system, then
# R_0 = (1)(\Theta_0) + (0)(\Theta_1} + (0)(\Theta_2)
# and
# Z_0 = R_0 - B_0
#
# This allows one to define restraints that are linear combinations
# of collective variables.  For cases in which only isolated collective variables
# are restrained, the matrix C_{k1} is the NxN identity matrix where N is
# the number of collective variables.
# 
proc cvc_getcvi { cvc } {
    set r {}
    set i 0
    foreach cv $cvc {
	if {[expr $cv > 0.0]} {
	    lappend r $i
	}
	incr i
    }
    return $r
}



proc rlist_setup { ncv } {
  global ds
  global restr
  global temperature
  global TAMDverbose

  for {set i 1} {$i<=$restr(num)} {incr i} {

    # Set cvc
    set restr($i.cvc) {}
    for {set j 0} { $j < $ncv} {incr j} {lappend restr($i.cvc) 0}
    foreach cvindex [lsort -unique $restr($i.cv)] {
      set j [expr $cvindex-1]
      set restr($i.cvc) [lreplace $restr($i.cvc) $j $j 1]
    }


    # Defaults and consistency in function
    if {![info exists restr($i.func)]} {set restr($i.func) HARMONIC}
    set restr($i.func) [string toupper $restr($i.func)]
    if {![info exists restr($i.k)]} {print "CFACV) WARNING: restraint with spring constant 0"; set restr($i.k) 0}


    # Defaults and consistency in boundaries
    if {![info exists restr($i.bound)]} {set restr($i.bound) NADA}
    set restr($i.bound) [string toupper $restr($i.bound)]
    switch $restr($i.bound) {
      NADA {
        if {![info exists restr($i.boundk)]} {set restr($i.boundk) 0.}
      }
      PERIODIC {
        set restr($i.min) [expr -1*acos(-1)]
        set restr($i.max) [expr acos(-1)]
        if {![info exists restr($i.boundk)]} {set restr($i.boundk) 0.}
      }
      PBC {
        if {![info exists restr($i.min)]}    {error "CFACV) ERROR: restraint need max value"}
        if {![info exists restr($i.max)]}    {error "CFACV) ERROR: restraint need min value"}
      }
      SOFTUPPER {
        if {![info exists restr($i.max)]}    {error "CFACV) ERROR: restraint need max value"}
        if {![info exists restr($i.boundk)]} {error "CFACV) ERROR: restraint need boundk value"}
      }
      SOFTLOWER {
        if {![info exists restr($i.min)]}    {error "CFACV) ERROR: restraint need min value"}
        if {![info exists restr($i.boundk)]} {error "CFACV) ERROR: restraint need boundk value"}
      }
      SOFT {
        if {![info exists restr($i.min)]}    {error "CFACV) ERROR: restraint need max value"}
        if {![info exists restr($i.max)]}    {error "CFACV) ERROR: restraint need min value"}
        if {![info exists restr($i.boundk)]} {error "CFACV) ERROR: restraint need boundk value"}
      }
      default {error "ERROR: $restr($i.bound) is not a valid boundary type."} 
    }

    if {![info exists restr($i.boundk)]} {set restr($i.boundk) 0.}
    if {![info exists restr($i.min)]}    {set restr($i.min) 0.}
    if {![info exists restr($i.max)]}    {set restr($i.max) 0.}
        
   
    # Output
    if {[info exists restr($i.outfile)]} {
      if {![info exists restr($i.outfreq)]}  {error "CFACV) ERROR: restraint need output frequency"}
    } else {
      if {[info exists restr($i.outfile)]}  {error "CFACV) ERROR: restraint need output file"}
      set restr($i.outfile) ""
      set restr($i.outfreq) -1
    }

    # Allocating ds->restr[i] in the C code.
    set pcvc [ListToArray $restr($i.cvc)]
    set inicial 0 ; # obsolet when r->z=r->val is placed at first iteration (see cfacv.c)
    set restr($i.address) [DataSpace_AddRestr $ds $restr($i.k) $inicial $ncv $pcvc $restr($i.func) $restr($i.min) $restr($i.max) $restr($i.bound) $restr($i.boundk) $restr($i.outfile) $restr($i.outfreq)]
     
    # Default and consistency in type. Allocating ds->restr[i]->tamdOpt in the C code.
    if {![info exists restr($i.type)]} {error "CFACV) ERROR: restraint need a type"}
    set restr($i.type) [string toupper $restr($i.type)]
    switch $restr($i.type) {
      TAMD {
        if {![info exists restr($i.g)]}     {error "CFACV) ERROR: TAMD restraint need g"}
        if {![info exists restr($i.dt)]}    {error "CFACV) ERROR: TAMD restraint need dt"}
        if {![info exists restr($i.chid)]}  {set restr($i.chid) 0}
        if {![info exists restr($i.chdm)]}  {set restr($i.chdm) 0}

        if {[info exists restr($i.temp)]}  {
          set ktaux [expr $restr($i.temp)*0.001987191]
        }

        if {[info exists restr($i.biasf)]}  {

          if {![info exists temperature]} {error "CFACV) ERROR: restr(biasf) declared but not temperature"}
          set ktaux [expr $temperature*$restr($i.biasf)*0.001987191]

          if {[info exists restr($i.temp)]} {
            if [expr {abs($temperature-$restr($i.temp)) < 1e-6}] {error "CFACV) ERROR: Declared restr(temp) differ from the desire restr(biasf)."}
          }
        }

        if {![info exists ktaux]}  {error "CFACV) ERROR: TAMD temperature info is needed."}

        # DataSpace_AddTamdOpt $ds $ir $restr($i.g) [expr $temperature*$restr($i.biasf)*0.001987191] $restr($i.dt)
        restr_AddTamdOpt $restr($i.address) $restr($i.g) $ktaux $restr($i.dt) $restr($i.chid) $restr($i.chdm)
      }
      SMD {
        if {![info exists restr($i.target)]} {error "CFACV) ERROR: SMD restraint need target"}
        if {![info exists restr($i.ini)]} {error "CFACV) ERROR: SMD restraint need ini"}
        if {![info exists restr($i.end)]} {error "CFACV) ERROR: SMD restraint need end"}
        restr_AddSmdOpt $restr($i.address) $restr($i.target) $restr($i.ini) $restr($i.end)
      }
      default {error "ERROR: $restr($i.type) is not a valid boundary type."} 
    } 

  }

  if {[info exists TAMDverbose]} {foreach ri [array get restr] {print "CFACV) $ri"}}
  
}



proc my_getcellsize { XSCFILE } {
    if {[info exists XSCFILE] && [file exists $XSCFILE]} {
    } else {
      error "CFACV) ERROR: no xsc file for box info."
    }
    #print "CFACV) DEBUG my_getcellsize opening $XSCFILE"
    set fp [open $XSCFILE "r"]
    set dat [read -nonewline $fp]
    #print "CFACV) DEBUG my_getcellsize dat:"
    #print "$dat"
    set lines [split $dat \n]
    #print "CFACV) DEBUG my_getcellsize lines:"
    #print "$lines"
    set lastline [lindex $lines end]
    lappend L [lindex $lastline 1]
    lappend L [lindex $lastline 5]
    lappend L [lindex $lastline 9]
    print "CFACV) INFO Box size $L"
    close $fp
    return $L
}

proc my_getorigin { XSCFILE } {
    if {[info exists XSCFILE] && [file exists $XSCFILE]} {
    } else {
      error "CFACV) ERROR: no xsc file for box info."
    }
    set fp [open $XSCFILE "r"]
    set dat [read -nonewline $fp]
    set lines [split $dat \n]
    set lastline [lindex $lines end]
    lappend O [lindex $lastline 10]
    lappend O [lindex $lastline 11]
    lappend O [lindex $lastline 12]
    print "CFACV) INFO origin $O"
    close $fp
    return $O
}

proc initialize_base_template { pdb } {
    mol new  $pdb

    set molid [molinfo top get id]

    set all [atomselect top all]
    $all set beta 0
    $all set occupancy 1.0

    $all delete
   
    return $molid
}

proc within_range {selList selB r} {
  set nc 0
  foreach sel $selList {
      set ABlist [measure contacts $r $sel $selB]
      set Alist [lindex $ABlist 0]
      set Blist [lindex $ABlist 1]
      if {[llength $Alist] > 0} {
          incr nc [llength $Alist]
      }
  }
  return $nc
}
     
