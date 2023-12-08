# tcl library for CFACV
# Cameron F Abrams 2009-14
# Sergio A Paz 2014-15

# This script is intended to be sourced before namd tclforce feature.  Therefore
# some of the procedures or variables used here are intrinsec of NAMD or can be
# found in the NAMD config file or in NAMD tcl master script for replica
# exechange.

# load the C-module
load ${CFACV_BASEDIR}/lib/libcfacv.so cfa_cvlibc

# load tcl procedures
source ${CFACV_BASEDIR}/include/data.tcl
source ${CFACV_BASEDIR}/include/restraints.tcl
source ${CFACV_BASEDIR}/include/cvs.tcl
source ${CFACV_BASEDIR}/include/chapeau.tcl

# BEGIN FUNCTION DEFINITION ###############################

proc cfacv_banner { argv } {
    global CFACV_VERSION
    global CFACV_BASEDIR
    print "================================================"
    print "CFACV) Collective Variables in VMD/NAMD v. $CFACV_VERSION"
    cfacvBanner
    print "CFACV) 2009-2014, Cameron F Abrams, Drexel University"
    print "CFACV) 2014-2015, Sergio A Paz, Drexel University"
    print "CFACV) 2016-2019, Sergio A Paz, Universidad Nacional de Córdoba - INFIQC/CONICET"
    print "================================================"
    print "Base directory: $CFACV_BASEDIR"
    print "argv: $argv"
    print "================================================"
    flush stdout
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


proc Tcl_UpdateDataSpace { ds lC groups first timestep bias } {
    upvar $lC p

    # Move group center positions to dataspace
    set i 0
    foreach g $groups {
      ListToArray_Data [DataSpace_centerPos $ds $i] $p($g)
      incr i
    }
    MyParanoiaCheck $ds "tripped after moving data to dataspace"

    # Compute CV's within dataspace
    set localbias [DataSpace_ComputeCVs $ds]
    set bias [expr $localbias+$bias]
    MyParanoiaCheck $ds "tripped after computing CV's"
    # At this point, the (x,y,z) position data is no longer needed.
    # We can now write into its space the (x,y,z) restraint forces.

    # Compute restraining forces
    DataSpace_RestrainingForces $ds $first $timestep $bias   
    MyParanoiaCheck $ds "tripped after computing restraining forces"

    # Transmit forces back to groups
    set i 0
    foreach g $groups {
        addforce $g [ArrayToList [DataSpace_centerPos $ds $i] 3]
        incr i
    }
    # puts "FUERZAS g1 [ArrayToList [DataSpace_centerPos $ds 0] 3]" 

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


# END FUNCTION DEFINITION ####################



# Global variables
set CFACV_VERSION 0.30
set kB_kcm 0.001987191
set first 1
set masses {}

# Banner
cfacv_banner NAMD

# Set the pseudorandom number generator seed
if {![info exists seed]} {
  print "CFACV) setting seed from clock"
  set seed [clock clicks]
}
print "CFACV) random seed $seed"

# Read the template PDB file that identifies subdomain memberships
parse_pdb $labelPDB pdb
print "CFACV) nCenters $pdb(nctrs)"

# Allocate dataspace in C code
set ds [NewDataSpace $pdb(nctrs) $cv(num) $restr(num) $seed]
 
# Add CVs to the datapsace structure
cvs_setup pdb
 
# Add restraints to the dataspace structure
rlist_setup $cv(num)

# Add chapeau to the dataspace structure 
# TODO: improve chapeau global variable
if {[info exists CFACV_doAnalyticalCalc]} {
  if {$CFACV_doAnalyticalCalc == 1} {
    if {![info exists USETAMDFORCES]} {set USETAMDFORCES 0}

    if {$restr(num) == 0} {error "No restraints to set a dataspace chapeau"}
    
    # Replica exechange mode: 
    if {![info exists NUMREP]} {set NUMREP 1}
    # TODO: NUMREP podria ser un numero de chapeaus distintos (estilo el
    # problema sodio cloro) a obtener en una simulacion. En ese caso
    # corregir.

    # Default modes
    if {![info exists DIMEN]} {set DIMEN 1}
    if {![info exists PERIODIC]} {
      for {set i 0} {$i < $DIMEN} {incr i} {lappend PERIODIC 0}
    }

    # Saving for allocate chapeau functions
    chapeau_setup $NUMREP $ds $DIMEN $SPLINEMIN $NKNOTS $CUTOFF $PERIODIC $BEGINEVOLVEPARAMETERS  $BEGINSOLVELAM $USETAMDFORCES $BINREPORTPARAMFILE $BINREPORTPARAMFREQ $BINOUTPUTLEVEL $LAMUPDATEINTERVAL
  } else {
    if {[info exists CUTOFF]} {error "Chapeau are not initialized"}
    if {[info exists DIMEN]} {error "Chapeau are not initialized"}
    if {[info exists SPLINEMIN]} {error "Chapeau are not initialized"}
    if {[info exists NKNOTS]} {error "Chapeau are not initialized"}
    if {[info exists PERIODIC]} {error "Chapeau are not initialized"}
  }
} else {
  if {[info exists CUTOFF]} {error "Chapeau are not initialized"}
  if {[info exists DIMEN]} {error "Chapeau are not initialized"}
  if {[info exists SPLINEMIN]} {error "Chapeau are not initialized"}
  if {[info exists NKNOTS]} {error "Chapeau are not initialized"}
  if {[info exists PERIODIC]} {error "Chapeau are not initialized"}
} 

# Some default values regarding output
if {![info exists TAMDof]}          {set TAMDof 1 }
if {![info exists TAMDbinof]}       {set TAMDbinof 1 }
if {![info exists TAMDoutputlevel]} {set TAMDoutputlevel 3}

set TAMDbinOutputFileFP 0
if {[info exists TAMDbinOutputFile]} {
    set TAMDbinOutputFileFP [my_binfopen $TAMDbinOutputFile "w" $TAMDoutputlevel $ds]
    print"CFACV) Binary TAMD output to $TAMDbinOutputFile"
}



