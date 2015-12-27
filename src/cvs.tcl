
proc read_cvs { cv_file cv_list pdb} {
  #and then cv_list is the append of dat.... why all the stuff?
  upvar $pdb pdbline 
  upvar $cv_list cvL

  set cvL {}

  # read the cv file
  set inStream [open $cv_file "r"]
  set data [read -nonewline $inStream]
  close $inStream
  set lines [split $data \n]

  set ncv 0
  foreach line $lines {

    # Allow empty lines and commentries
    if {[string index $line 0] == "#"} continue
    if {[string length $line] < 0} continue

    # Supouse that the line is "CARTESIAN_X   0"...
    set dat {}
    foreach l $line { lappend dat $l } ;#... dat is "CARTESIAN_X 0"
    set typ [lindex $dat 0]            ;#... typ is "CARTESIAN_X"
    set ind [lreplace $dat 0 0]        ;#... ind is "0"

    # Search for optional strings after the index list
    set optlist {}
    for {set i 0} {[string is integer [lindex $ind $i]] && $i < [llength $ind]} { incr i } {}
    if {$i < [llength $ind]} {
        # Put in optlist that string and let in ind only the numbers
        set optlist [lreplace $ind 0 [expr $i-1]]
        set ind [lreplace $ind $i end]
    }
    set nInd [llength $ind]
    set nOpt [llength $optlist]
    
    # Si no hay indices, busco palabras clave de seleccion
    if {$nInd==0} {
      set j 0
      foreach o $optlist {
        incr j

        # fill index by resname
        if {$o=="resname"} {
          print [lindex $optlist  $j]
          set resn [lindex $optlist  $j]
          foreach l $pdbline {
            set index [lindex $l 4]
            if {$resn==[lindex $l 2]} {lappend ind [expr $index-1]}
          }

          # Remove resname from option list
          set j [expr $j-1]
          set optlist [lreplace $optlist $j $j ]
          set optlist [lreplace $optlist $j $j ]

          break
        }

      }
    }
    set nInd [llength $ind]

    # Si todavia no hay indices, doy error
    if {$nInd==0} {error "Error, no indices para la linea $line"}


    lappend cvL [list $typ $ind $optlist] ;#... cvL is cvL+dat
    incr ncv
    
    # Check for a valid CV
    switch $typ {
      ZSDCIRCLE {
        set_zsd_circle [lindex $optlist 0] [lindex $optlist 1]  [lindex $optlist 2] [lindex $optlist 3]
      }
      ZSDXRANGE {
        set_zsd_circle [lindex $optlist 0] 0.0  [lindex $optlist 1] [lindex $optlist 2]
      }
      ZSDRING {
        set_zsd_ring [lindex $optlist 0] [lindex $optlist 1]  [lindex $optlist 2] [lindex $optlist 3] [lindex $optlist 4] 
      }
      BOND -
      S -
      ANGLE -
      DIHED -
      COGX -
      COGY -
      COGZ -
      CARTESIAN_X -
      CARTESIAN_Y -
      CARTESIAN_Z {}
      default {error "ERROR: $typ is not a valid CV type in $cv_file."}
    }
  }

  return $ncv
}

proc generate_cv { cvOpt p } {

    puts -nonewline "CFACV) Generating cv.inp from $cvOpt..."; flush stdout
    set retList {}
    set segidList {}

    set nCntr [llength $p]
    
    foreach s $p {
	set segid [lsort -unique [$s get segid]]
	if {[llength $segid] > 1} {
	    print "ERROR: a center owns atoms from more than one segid: $segid"
	}
	lappend segidList $segid
    }

#    print "DB: generate_cv: segidList: $segidList" 

    if {($cvOpt == "autogenerate-bonds" || $cvOpt == "all-to-all-bonds")} {
	for {set ii 0} {$ii < [expr $nCntr-1]} { incr ii } {
	    for {set jj [expr $ii + 1]} { $jj < $nCntr } { incr jj } {
		set thisind [list $ii $jj]
		lappend retList [list "BOND" $thisind]
	    }
	}
    } elseif {[string first "shortest" $cvOpt]!=-1} {
	set tl [split $cvOpt -]
	set nIntraSeg [lindex $tl 1]
	if {[llength $tl] > 2} {
	    set nInterSeg [lindex $tl 2]
	}
	set ii 0
	set jj 0
	set intraSegL {}
	set interSegL {}
	foreach s1 $p {
	    set iseg [lindex $segidList $ii]
	    set c1 [measure center $s1 weight mass]
	    set jj [expr $ii + 1]
	    foreach s2 [lreplace $p 0 $ii] {
		set jseg [lindex $segidList $jj]
		set c2 [measure center $s2 weight mass]
		set d12 [veclength [vecsub $c1 $c2]]
		if {$iseg == $jseg} {
		    lappend intraSegL [list $ii $jj $d12]
		} else {
		    lappend interSegL [list $ii $jj $d12]
		}
		incr jj
	    }
	    incr ii
	}
	set intraSegL [lsort -real -index 2 $intraSegL]
	set interSegL [lsort -real -index 2 $interSegL]
	for {set i 0} {$i < $nIntraSeg} {incr i} {
	    set tl [lindex $intraSegL $i]
	    set ind [list [lindex $tl 0] [lindex $tl 1]]
	    lappend retList [list "BOND" $ind]
	}
	if {[info exists nInterSeg]} {
	    for {set i 0} {$i < $nInterSeg} {incr i} {
		set tl [lindex $interSegL $i]
		set ind [list [lindex $tl 0] [lindex $tl 1]]
		lappend retList [list "BOND" $ind]
	    }
	}
	set retList [lsort -index 1 $retList]
    } elseif {[string first "maxlength" $cvOpt] != -1} {
	set tl [split $cvOpt =]
	set maxL_intra [lindex $tl 1]
	if {[llength $tl] > 2} {
	    set maxL_inter [lindex $tl 2]
	}
	set ii 0
	set jj 0
	set intraSegL {}
	set interSegL {}
	foreach s1 $p {
	    set iseg [lindex $segidList $ii]
	    set c1 [measure center $s1 weight mass]
	    set jj [expr $ii + 1]
	    foreach s2 [lreplace $p 0 $ii] {
		set jseg [lindex $segidList $jj]
		set c2 [measure center $s2 weight mass]
		set d12 [veclength [vecsub $c1 $c2]]
		if {$iseg == $jseg} {
		    if {[expr $d12 < $maxL_intra]} {
			lappend intraSegL [list $ii $jj $d12]
		    }
		} else {
		    if {[expr $d12 < $maxL_inter]} {
			lappend interSegL [list $ii $jj $d12]
		    }
		}
		incr jj
	    }
	    incr ii
	}
	set intraSegL [lsort -real -index 2 $intraSegL]
	set interSegL [lsort -real -index 2 $interSegL]
	set nIntraSeg [llength $intraSegL]
	set nInterSeg [llength $interSegL]
	for {set i 0} {$i < $nIntraSeg} {incr i} {
	    set tl [lindex $intraSegL $i]
	    set ind [list [lindex $tl 0] [lindex $tl 1]]
	    lappend retList [list "BOND" $ind]
	}
	if {[info exists nInterSeg]} {
	    for {set i 0} {$i < $nInterSeg} {incr i} {
		set tl [lindex $interSegL $i]
		set ind [list [lindex $tl 0] [lindex $tl 1]]
		lappend retList [list "BOND" $ind]
	    }
	}
	set retList [lsort -index 1 $retList]
    } elseif {[string first "cartesian" $cvOpt] != -1} {
	for {set ii 0} {$ii < $nCntr} { incr ii } {
	    foreach c {X Y Z} {
		lappend retList [list "CARTESIAN_${c}" $ii]
	    }
	}
    } elseif {[string first "peptide_tilt" $cvOpt] != -1} {
	set nDihed [expr ($nCntr/2)-1]
	# O 1, 3, 5, 7, 9, etc.
        # H 2, 4, 6, 8, 10, etc.
	# dihed: 2-1-3-4, 4-3-5-6, 6-5-7-8, 8-7-9-12, etc.
	# but begin counting at 0!!!
	set p1 1
	set p2 0
	set p3 2
	set p4 3
	for {set i 0} {$i < $nDihed} {incr i} {
	    lappend retList [list "DIHED" $p1 $p2 $p3 $p4]
	    incr p1 2
	    incr p2 2
	    incr p3 2
	    incr p4 2
	}
    } elseif {[string first "phi-psi" $cvOpt] != -1} {
	set nPsi [expr ($nCntr/3)-1]
	print "CFACV) [expr ($nCntr/3)] residues -> $nPsi psi angles"
	set psia 0
	set psib 1
	set psic 2
	set psid 3
	set phia 2
	set phib 3
	set phic 4
	set phid 5
	for {set i 0} {$i < $nPsi} {incr i} {
	    lappend retList [list "DIHED" $psia $psib $psic $psid]
	    lappend retList [list "DIHED" $phia $phib $phic $phid]
	    incr psia 3
	    incr psib 3
	    incr psic 3
	    incr psid 3
	    incr phia 3
	    incr phib 3
	    incr phic 3
	    incr phid 3
	}
    } else {
	print "ERROR: cvOpt $cvOpt is not recognized."
    }
    print "done."
    return $retList
}

proc output_cv { cvList fileName } {
    set fp [open "$fileName" "w"]
    foreach cv $cvList {
	set typ [lindex $cv 0]
	set ind [lreplace $cv 0 0]
	print $fp "$typ $ind"
    }
    close $fp
}


proc cvs_setup { } {
  global cvList
  global ds

  foreach cv $cvList {
    set typ [lindex $cv 0]
    set pind [intListToArray [lindex $cv 1]]
    set nind [llength [lindex $cv 1]]

    # Assuming cv=CARTESIAN_X 0 then typ=CARTESIAN_X, pind={0}, nind=1
    DataSpace_AddCV $ds $typ $nind $pind 
  }

}
  
