# tcl library for CFACV
# Cameron F Abrams
# 2009-14
#
# Any Tcl script that uses this library must set the value
# of CFACV_BASEDIR correctly.
#
set CFACV_VERSION 0.30
# load the C-module
load ${CFACV_BASEDIR}/cfacv.so cfa_cvlibc

proc cfacv_banner { argv } {
    global CFACV_VERSION
    global CFACV_BASEDIR
    print "================================================"
    print "CFACV) Collective Variables in VMD/NAMD v. $CFACV_VERSION"
    cfacvBanner
    print "CFACV) 2009-2011, Cameron F Abrams, Drexel University"
    print "================================================"
    print "Base directory: $CFACV_BASEDIR"
    print "argv: $argv"
    print "================================================"
    flush stdout
}

###########################################################
# Primitive data handling procedures
###########################################################

# ListToArray: allocates a new array and assigns it values
# from the Tcl list $l; returns pointer to new array.  List
# elements are treated as double-precision by default.
#
proc ListToArray {l} {
    set length [llength $l]
    set a [new_array $length]
    set i 0
    foreach item $l {
        array_setitem $a $i $item
        incr i 1
    }
    return $a
}

# ListToArray_Data:  Assigns elements of an existing 
# array pointed to by $a from the tcl list $l
#
proc ListToArray_Data { a l } {
    set i 0
    foreach item $l {
	array_setitem $a $i $item
	incr i 1
    }
    return $a
}

# intListToArray: like listToArray except elements of
# $l are treated as integers
#
proc intListToArray {l} {
   set length [llength $l]
   set a [new_arrayint $length]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

# intListToArray_Data: list listToArray_Data except elements
# of $l are treated as integers
#
proc intListToArray_Data {a l} {
   set length [llength $l]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

# ArrayToList:  Creates a new tcl list of length $n and
# assigns it elements of existing array $a
#
proc ArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [array_getitem $a $i]
    }
    return $l
}

# intArrayToList: like ArrayToList but treats
# elements of $a as integers
#
proc intArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [arrayint_getitem $a $i]
    }
    return $l
}

proc getArg { argv key altkey def } {
    set i 0
    while {$i < [llength $argv] && [lindex $argv $i] != "-$key" && [lindex $argv $i] != "-$altkey"} {
#	print "DB: getArg considers [lindex $argv $i]"
	incr i
    }
    if {$i < [llength $argv] && ([lindex $argv $i] == "-$key" || [lindex $argv $i] == "-$altkey")} {
	incr i
	set rv {}
	for {} {$i < [llength $argv] && [string index [lindex $argv $i] 0] != "-"} {incr i} {
	    lappend rv [lindex $argv $i]
#	    print "DB: getArg building rv as $rv"
	}
	return $rv
    }
    return $def
}

###############################################################
# read_centersPDB
###############################################################
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
#  masses: a tcl list in which the i'th element
#          is the total mass of center i
#
#  return value is the number of centers
#
###############################################################
proc read_centersPDB { templatePdb serArr mass } {
    upvar $serArr serArray
    upvar $mass masses

    set inStream [open $templatePdb r]
    set data [read $inStream]
    close $inStream
    set lines [split $data \n]

    set n 0
    set indices {}
    foreach line $lines {
	set cardstr [string range $line 0 3]
	set occustr [string trim [string range $line 54 59]]
	set betastr [string trim [string range $line 60 65]]
	#print "o|$occustr|b|$betastr|"
	if {([string equal $cardstr "ATOM"]) && [expr 0==[string equal $betastr "0.00"]]} {
	    set t [expr int($betastr)-1]
#	    print "DB: found betastr |$betastr| => group $t"
	    set index [lsearch $indices $t]
	    if { [expr $index == -1] } {
		lappend indices $t
	    }
	}
	incr n
    }
    set indices [lsort -integer $indices]
    set nMon [llength $indices]

    print "DB: read_centersPDB: indices $indices nMon $nMon"

    for {set i 0} {$i < $nMon} {incr i} {
	lappend serArray {}
	lappend masses 0.0
    }
    set n 0
    foreach line $lines {
	set cardstr [string range $line 0 3]
	set serlstr [string trim [string range $line 5 11]]
        # occustr is used to hold atomic mass in amu.
	set occustr [string trim [string range $line 54 59]]
	set betastr [string trim [string range $line 60 65]]
	# if the beta field is nonzero, it is assumed to be an integer
	# that indicates which center this atom belongs to
	if {([string equal $cardstr "ATOM"]) && [expr 0==[string equal $betastr "0.00"]]} {
	    # identify which center to which this atom belongs
	    set t [expr int($betastr)-1]
	    # extract that center's current list of atom indices from the list of lists
	    set nlist [lindex $serArray $t]
	    # add this index to the list; index is serial-1

	    # since begining with 100000, serial numbers are stored and
	    # read in hex, need to check whether this is a hex number
	    # -- I'll do this by seeing if it is not an integer
#	    print "adding serl $serlstr to $nlist"
	    if {[string is integer $serlstr]} {
		set serl [expr $serlstr]
	    } else {
		set serl [expr 0x$serlstr]
	    }
	
	    lappend nlist $serl
	    # return the newly appended list of indices to the list of lists
	    set serArray [lreplace $serArray $t $t $nlist]
	    # update the mass of this pseudoatom
	    set am [expr double($occustr)]
	#    print "DB: read_centersPDB: atom serial $serlstr mass $am to mon $t"
	    set masses [lreplace $masses $t $t [expr [lindex $masses $t] + $am]]
	}
	incr n
    }
    print "DB: returning nMon $nMon"
#    print "DB: serArray : $serArray"
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

######################################################################
# read_centersVMD
#
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
#
#####################################################################
proc read_centersVMD { pL mL fileName molID } {
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
#	print "DB: making selection from residues $il"
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

proc read_cvs { cv_file cv_list nMon } {
    upvar $cv_list cvL

    set cvL {}

    set inStream [open $cv_file "r"]
    set data [read -nonewline $inStream]
    close $inStream
    set lines [split $data \n]
    set ncv 0
    foreach line $lines {
	if {[string index $line 0] != "#" && [string length $line] > 0} {
	    set dat {}
	    foreach l $line { lappend dat $l }
	    set typ [lindex $dat 0]
	    set ind [lreplace $dat 0 0]
	    set optlist {}
	    for {set i 0} {[string is integer [lindex $ind $i]] && $i < [llength $ind]} { incr i } {}
	    if {$i < [llength $ind]} {
		set optlist [lreplace $ind 0 [expr $i-1]]
		set ind [lreplace $ind $i end]
	    }
	    set nInd [llength $ind]
	    set nOpt [llength $optlist]
	    
#	    print "DB: $nInd: ind     $ind"
#	    print "DB: $nOpt: optlist $optlist"
	    
	    lappend cvL [list $typ $ind $optlist]
	    incr ncv
	    
	    switch $typ {
		BOND -
		ANGLE -
		DIHED -
		CARTESIAN_X -
		CARTESIAN_Y -
		CARTESIAN_Z {
		}
		default {
		    print "ERROR: $typ is not a valid CV type in $cv_file."
		    exit
		}
	    }
	}
    }

#    foreach cv $cvL {
#	print "DB: cv $cv"
#	print "DB: typ: [cv_gettyp $cv]"
#	print "DB: ind: [cv_getind $cv]"
#	print "DB: optlist: [cv_getoptlist $cv]"
#    }

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

####################################################################
# Procedures for handling restraints
####################################################################
#
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

proc restr_getcvc { r } {
    return [lindex $r 0]
}
proc restr_getoptlist { r } {
    return [lindex $r 1]
}

proc restr_getopt { r key altkey def } {
    set optlist [restr_getoptlist $r]
    foreach opt $optlist {
#	print "DB: querying opt $opt for $key..."
	set thiskey [lindex $opt 0]
	if {$thiskey == $key || $thiskey == $altkey} {
	    return [lindex $opt 1]
	}
    }
    return $def
}

proc read_restraints { restr_file ncv restrList } {
    upvar $restrList rL

    set inStream [open $restr_file "r"]
    set data [read -nonewline $inStream]
    close $inStream
    set lines [split $data \n]

    set rL {}
    set nR 0
    foreach line $lines {
	if {[string index $line 0] != "#" && [string length $line] > 0} {
	    set cvc [lreplace $line $ncv end]
	    set optlist [lreplace $line 0 [expr $ncv - 1]]
	   # print "DB: cvc $cvc"
	    #print "DB: optlist $optlist"
	    lappend rL [list $cvc $optlist]
	    incr nR
	}
    }
    #print "DB: rL $rL"
    return $nR
}

proc create_single_cv_restraints { ncv restrList restrPARAMS } {
    upvar $restrList rL
    set rL {}
    set nR 0
#    set cvc [lrepeat $ncv 0]
    set cvc {}
    for {set i 0} { $i < $ncv} {incr i} {
	lappend cvc 0
    }
    for {set i 0} { $i < $ncv} {incr i}  {
	set tcvc [lreplace $cvc $i $i 1]
	lappend rL [list $tcvc $restrPARAMS]
    }
    return $ncv
}

proc Tcl_NewDataSpace { nC cvL rL seed } {
    global doPairCalc
    set nCV [llength $cvL]
    set nR  [llength $rL]

    set ds [NewDataSpace $nC $nCV $nR $seed]

    foreach cv $cvL {
	set typ [lindex $cv 0]
	set pind [intListToArray [lindex $cv 1]]
	set nind [llength [lindex $cv 1]]
	DataSpace_AddCV $ds $typ $nind $pind 
    }

    foreach r $rL {
	set cvc [restr_getcvc $r]
	set nCV [llength $cvc]
	set pcvc [ListToArray $cvc]
	set k    [restr_getopt $r "SpringConstant" "k" 0.00]
	set targ [restr_getopt $r "TargetValue"    "b" 0.00]
	set RestraintFunction [string toupper [restr_getopt $r "RestraintFunction" "rf"  HARMONIC]]
	set zmin [restr_getopt $r "Minimum" "min" 0.00]
	set zmax [restr_getopt $r "Maximum" "max" 0.00]
	if { $RestraintFunction == "PERIODIC" } {
	    set zmin [expr -1*acos(-1)]
	    set zmax [expr acos(-1)]
	}
	set ir [DataSpace_AddRestr $ds $k $targ $nCV $pcvc $RestraintFunction $zmin $zmax]
	set TAMDgamma [restr_getopt $r "TAMDgamma" "gamma" -1.0]
	set TAMDkT    [restr_getopt $r "TAMDkT"    "kT"    -1.0]
	set TAMDdt    [restr_getopt $r "TAMDdt"    "dt"    -1.0]
	if {$TAMDgamma != -1 && $TAMDkT != -1 && $TAMDdt != -1} {
	    DataSpace_AddTamdOpt $ds $ir $TAMDgamma $TAMDkT $TAMDdt
	}
	set SMDTarget     [restr_getopt $r "SMDTargetValue" "smd_z"     -1.0]
	set SMDInitStep   [restr_getopt $r "SMDInitialStep" "smd_t0"    -1.0]
	set SMDFinalStep  [restr_getopt $r "SMDFinalStep"   "smd_t1"    -1.0]
	if {$SMDTarget != -1 && $SMDInitStep != -1 && $SMDFinalStep != -1} {
	    DataSpace_AddSmdOpt $ds $ir $SMDTarget $SMDInitStep $SMDFinalStep
	}
    }

    return $ds
}

proc my_getcellsize { XSCFILE } {
    if {[info exists XSCFILE] && [file exists $XSCFILE]} {
    } else {
      print "CFACV) ERROR: no xsc file for box info."
      exit
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
      print "CFACV) ERROR: no xsc file for box info."
      exit
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

proc Tcl_InitializePairCalc { ds XSCFILE cutoff nlcutoff begin_evolve usetamdforces reportparamfreq spline_min nKnots splineoutputfile splineoutputfreq splineoutputlevel updateinterval } {
    set LL [my_getcellsize $XSCFILE]
    set O [my_getorigin $XSCFILE]
    print "CFACV) DEBUG: Tcl_InitializePairCalc box size [lindex $LL 0] [lindex $LL 1] [lindex $LL 2]"
    print "CFACV) DEBUG: Tcl_InitializePairCalc origin   [lindex $O 0] [lindex $O 1] [lindex $O 2]"
    DataSpace_SetupPairCalc $ds [lindex $O 0] [lindex $O 1] [lindex $O 2] [lindex $LL 0] [lindex $LL 1] [lindex $LL 2] $cutoff $nlcutoff $begin_evolve $usetamdforces $reportparamfreq $spline_min $nKnots $splineoutputfile $splineoutputfreq $splineoutputlevel $updateinterval
    print "CFACV) INFO: Pair calc initialized"
}

proc Tcl_Reinitialize { ds restartINP } {
    if {[file exists $restartINP]} {
	print "CFACV) Reinitializing Z values from $restartINP"
	set fp [open $restartINP "r"]
	set zinp [read -nonewline $fp]
	close $fp
	DataSpace_SetRestraints $ds [ListToArray $zinp]
	return 0
    } else {
	print "CFACV) Error: could not find restart file $restartINP"
	print "CFACV) Initial z-values taken from CV's."
	return 1
    }
}

proc Tcl_UpdateDataSpace { ds lC groups first timestep } {
    upvar $lC p

    # Move group center positions to dataspace
    set i 0
    foreach g $groups {
	ListToArray_Data [DataSpace_centerPos $ds $i] $p($g)
	incr i
    }
    MyParanoiaCheck $ds "tripped after moving data to dataspace"

    # Compute CV's within dataspace
    DataSpace_ComputeCVs $ds
    MyParanoiaCheck $ds "tripped after computing CV's"
    # At this point, the (x,y,z) position data is no longer needed.
    # We can now write into its space the (x,y,z) restraint forces.

    # Compute restraining forces
    DataSpace_RestrainingForces $ds $first $timestep    
    MyParanoiaCheck $ds "tripped after computing restraining forces"

    # Transmit forces back to groups
    set i 0
    foreach g $groups {
	addforce $g [ArrayToList [DataSpace_centerPos $ds $i] 3]
	incr i
    }

    # Add restraint energy to NAMD energy structure
    addenergy [DataSpace_RestraintEnergy $ds]
}

proc Tcl_ObserveDataSpace { ds cntrSel frame } {
    set i 0
    foreach sel $cntrSel {
	$sel frame $frame
#	print "DB: Tcl_ObserveDataSpace center $i [measure center $sel weight mass]"
	ListToArray_Data [DataSpace_centerPos $ds $i] [measure center $sel weight mass]
	incr i
    }
    DataSpace_ComputeCVs $ds
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

# This routine uses the DataSpace_checkdata function
# to detect any "infs" or "nans" in the position data
# and dies with output if any are detected.
# Note this is a stub if a variable name PARANOIA
# is *not* defined in either the NAMD config file
# or the cfacv_tclforces script.
proc MyParanoiaCheck {ds msg} {
    global PARANOIA
    if {[info exists PARANOIA]} {
	if {$PARANOIA} {
	    if {[DataSpace_checkdata $ds]} {
		print "CFACV/PARANOIA) $msg"
		DataSpace_dump $ds
		exit
	    }
	}
    }
}

proc Tcl_DataSpace_InitKnots { ds filename } {
    if {[file exists $filename]} {
	DataSpace_InitKnots $ds $filename
    } {
	print "CFACV) error: initial knots file $filename not found."
	exit
    }

}