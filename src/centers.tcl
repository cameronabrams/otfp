
proc read_centersPDB { templatePdb serArr mass pdb} {
   
  # Sets up the centers according to the template PDB file
  # in which each atom's beta field designates the index of
  # the center to which it belongs; each atom can belong
  # to one and only one center!  Also, the occupancy
  # field should contain the atomic mass in amu.
  #
  # IMPORTANT: A beta of i means the atom belongs to the (i-1)th 
  # center! A beta of "0" indicates that the atom does not belong 
  # to a center.
  # 
  # returns:
  #  serArr: a tcl list in which the i'th element
  #          is a list of atom serial numbers (1-based indicies)
  #          belonging to center i
  #
  #          e.g.
  #          {1 2 4    } ;#group 1
  #          {3 5 7    } ;#group 2
  #          {6 9 12 11} ;#group 3
  #
  #  masses: a tcl list in which the i'th element
  #          is the total mass of center i
  #
  #          e.g.
  #          {1.004 15.99 1.004       } #group 1
  #          {1.004 15.99 1.004       } #group 2
  #          {1.004 15.99 1.004 1.004 } #group 3
  #
  #  return value is the number of centers
 
  upvar $serArr serArray
  upvar $mass masses
  upvar $pdb pdbline 

  # pdbline is a list with all the properties from the pdb file
  # this allow to access to the resname or segname later, which will
  # allow to use resname or segname keywords in the cv declaration

  # Open and read the pdb. Put all it in lines variable.
  set inStream [open $templatePdb r]
  set data [read $inStream]
  close $inStream
  set lines [split $data \n]

  # Get list indices with the center numbers and Nmon as the number of
  # centers
  set indices {}
  foreach line $lines {
      set cardstr [string range $line 0 3]
      set betastr [string trim [string range $line 60 65]]
      # if the beta field is nonzero, it is assumed to be an integer
      # that indicates which center this atom belongs to
      if {([string equal $cardstr "ATOM"]) && [expr 0==[string equal $betastr "0.00"]]} {
        lappend indices [expr int($betastr)-1]
      }
  }
  set indices [lsort -integer -unique $indices]
  set nMon [llength $indices]

  # contruct the void list of masses, atom serial and reslist
  for {set i 0} {$i < $nMon} {incr i} {
      lappend serArray {}
      lappend pdbline {}
      lappend masses 0.0
      lappend reslist 0
  }

  # append list of masses and serials ID to each group index
  set n 0
  foreach line $lines {
      set cardstr [string range $line 0 3]
      set serlstr [string trim [string range $line 5 11]]
      set resname [string trim [string range $line 17 20]]
      set occustr [string trim [string range $line 54 59]] ;# used to hold atomic mass in amu.
      set betastr [string trim [string range $line 60 65]]
      if {([string equal $cardstr "ATOM"]) && [expr 0==[string equal $betastr "0.00"]]} {
          
          set t [expr int($betastr)-1] ; #Center index
          set nlist [lindex $serArray $t] ; # get atom index list for center $t

          # add this index to the list; index is serial-1
          # since begining with 100000, serial numbers are stored and
          # read in hex, need to check whether this is a hex number
          # -- I'll do this by seeing if it is not an integer
#         print "adding serl $serlstr to $nlist"
          if {[string is integer $serlstr]} {
      	set serl [expr $serlstr]
          } else {
      	set serl [expr 0x$serlstr]
          }
          lappend nlist $serl
          # return the newly appended list of indices to the list of lists
          set serArray [lreplace $serArray $t $t $nlist]

          # update all the properties. Now, if is a group, pdbline
          # correspond to the last line of the group. TODO I don't know it
          # have sense for a group this mechanism, but it is good to
          # select many atoms in the cv.  I guess that it should be better
          # to put one entrance per atom in the pdbline, then the info of
          # the groups will be saved, and avoid serarray and any other
          # stuff that complicate the legibility
          set pdbline [lreplace $pdbline $t $t "$cardstr $serlstr $resname $occustr $betastr"]

          # update the mass of this pseudoatom
          set am [expr double($occustr)]
          set masses [lreplace $masses $t $t [expr [lindex $masses $t] + $am]]

          # Renames
          set reslist [lreplace $reslist $t $t $resname]

      }
      incr n
  }

  # FIXME: This code was here to allow compute chapeau functions separatedly
  # for different pair types of particles. For instance, this allow to
  # recover SOD SOD, CLA CLA and SOD CLA pair potentials in 1 TAMD
  # simulation. Each index has a number in ch_id which allow to sort the pair
  # in the different chapeau objects on the c code.  From the studies with
  # SOD CLA, this pair potentials will be OK only if the ficticius
  # temperature is the same that the real one.  On the other hand, a better
  # way to achive this is needed (without saving a lot of numbers in ch_id).
  # For understand how this worked, see the previous versions of the code.
  # set k 0
  # set test [expr 0==[string equal $chlist "{}"]]
  # for {set i 0} {$i < $nMon} {incr i} {
  #   for {set j [expr $i+1]} {$j < $nMon} {incr j} {

  #     if $test {
  #       set pair1 "[lindex $reslist $i] [lindex $reslist $j]" ;# e.g. "SOD CLA"
  #       set pair2 "[lindex $reslist $j] [lindex $reslist $i]" ;# e.g. "CLA SOD"
  #       set a [lsearch $chlist $pair1]
  #       set b [lsearch $chlist $pair2]
  #       
  #       set k [expr {$a>$b? $a: $b}]
  #     }

  #     # with the next expresion the [expr $j+($nCntr-2)*$i-($i-1)*$i/2-1] element of ch_id is
  #     # -1 for non intereset pair, or the corresponding n-index or the $chlist.
  #     lappend ch_id $k
  #   }
  # }

  # if {([string equal [lsort -unique $ch_id] "-1"])} {
  #   error "ERROR: No residues found to accelerate"
  # }

  print "DB: returning nMon $nMon"
  return $nMon
}


proc center_selections { molID serArray } {
    set p {}
    set n [llength $serArray]
    print "DB: center_selections thinks there are $n centers"
    set i 0
    for {set i 0} { $i < $n } { incr i } {
	set arr [lindex $serArray $i]
	set a [atomselect $molID "serial $arr"]
        $a global
	lappend p $a
    }
    return $p
}


proc read_centersVMD { pL mL fileName molID } {
  # If a molecule is loaded in VMD and a "centers.dat" file
  # (generated by mk_tPDB.tcl) exists, then it is possible
  # to define the atomselections for each center without
  # reading a template PDB.  Because "atomselections" are
  # not part of the Tcl interface to NAMD, this approach
  # is not useful for tcl_forces.  
  #
  # pL is a Tcl list of atomselections; this procedure only appends it.
  # fileName is the name of the "centers.dat" file
  # molID is the molecule id number of the molecule for
  # which centers are defined.
  #
  # The format of the "centers.dat" file is such that
  # one center is defined on each line.  The format
  # of a line is
  #
  # id mass residue_1 residue_2 ... residue_N
  #
  # where id is the index of the center, mass is its mass in amu
  # and residue_i is a string concatenation of the single-letter
  # amino acid designation and the residue number.
  upvar $pL p
  upvar $mL m
 
  set ofp [open "$fileName" "r"]
  set data [read -nonewline $ofp]
  close $ofp
  set lines [split $data \n]
  foreach l $lines {
      set il {}
      lappend m [lindex $l 1]
      set RL [lreplace $l 0 1]
      foreach ri $RL {
          lappend il [string range $ri 1 end]
      }
#     print "DB: making selection from residues $il"
      set a [atomselect $molID "protein and resid $il"]
      $a global
      lappend p $a 
  }
  return [llength $lines]
}

proc read_centers_residLists { ridL fileName } {
    upvar $ridL r
    set ofp [open "$fileName" "r"]
    set data [read -nonewline $ofp]
    close $ofp
    set lines [split $data \n]
    foreach l $lines {
	set il {}
	set RL [lreplace $l 0 1]
	foreach ri $RL {
	    lappend il [string range $ri 1 end]
	}
	lappend r $il 
    }
    return [llength $lines]    
}


proc assign_center_betas { p includeH } {
    set nCntr [llength $p]
    for {set i 0} { $i < $nCntr } { incr i } {
	set cntrIndex [expr $i + 1]
	set a [lindex $p $i]
	if { $includeH } {
	    $a set beta $cntrIndex
	} else {
	    set thisbeta {}
	    foreach e [$a get element] {
		if { $e != "H" } {
		    lappend thisbeta $cntrIndex
		} else {
		    lappend thisbeta 0.0
		}
	    }
	    $a set beta $thisbeta
	}
	$a set occupancy [$a get mass]
    }
}

proc report_centers { p fileName } {
    set nCntr [llength $p]
    set ofp [open "$fileName" "w"]
    for {set i 0} { $i < $nCntr } { incr i } {
	set a [lindex $p $i]
	set this_m [format "%8.2f" [vecsum [$a get mass]]]
	puts -nonewline $ofp "[format %3i $i] $this_m "
	foreach id [lsort -integer -unique -index 1 [$a get {resname resid chain}]] {
	    set indx "[aa_321 [lindex $id 0]][format %-3i [lindex $id 1]]"
	    puts -nonewline $ofp "$indx "
	}
	puts $ofp ""
    }
    close $ofp
}
              
 
