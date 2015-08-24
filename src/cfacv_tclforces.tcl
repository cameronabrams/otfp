# tclforces file for CFACV
# Cameron F Abrams
# 2009-14

# This script is intended to be sourced by namd tclforce feature.  Therefore
# some of the procedures or variables used here are intrinsec of NAMD or can be
# found in the NAMD config file or in NAMD tcl master script for replica
# exechange.
#
# Procedures or variables defined in the tcl master script for replica
# exechanges will be aviabales here. At the same time procedures or variables
# defined in the present scope will be aviabales in tcl master script only
# after first run command invocation are executed in that script.
# 
# NOTE: NAMD require that the tcl command "puts" be replaced by "print" here.
#

# Global variables
set CFACV_VERSION 0.30
set kB_kcm 0.001987191

# load the C-module
load ${CFACV_BASEDIR}/cfacv.so cfa_cvlibc

# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/cfacv.tcl

cfacv_banner NAMD

# Search for missing key
set tripped 0
foreach key {labelPDB cvINP} {
  if {![info exists $key]} {error "CFACV) ERROR: you must set $key in the NAMD config. file."}
}


#### Get the groups using the addgroup of tclforces

set serArray {}; # must have for addgroup
set pdbline {};
set masses {}

# Read the template PDB file that identifies subdomain memberships
set nCntr [read_centersPDB $labelPDB serArray masses pdbline]
#nCntr is the number of centers
#serArray is the atom serial list of each group
#masses is the atom mass list of each group
#ch_id holds the chapeau object for the ij pair type
print "CFACV) nCenters $nCntr  masses $masses"

# Set up the centers as "groups" for tclforces. This give 
# the list groups with the groups ids
set groups {}
for {set i 0} { $i < $nCntr } { incr i } {
    if {[info exists TAMDverbose]} {
	print "addgroup $i :"
	print "   [lindex $serArray $i]"
    }
    # addgroup is a NAMD command that return a "gN" id with N a small
    # integer
    set gid [addgroup [lindex $serArray $i]]
    lappend groups $gid
}


# Set the pseudorandom number generator seed
if {![info exists seed]} {
  print "CFACV) setting seed from clock"
  set seed [clock clicks]
}
print "CFACV) random seed $seed"

 
# Set up list of CV's
#cvList is some like {CARTESIAN_X 0} {CARTESIAN_Y 0} ....
set cvList {}
set nCV [read_cvs $cvINP cvList pdbline]
if {!$nCV} {error "CFACV) ERROR: wrong cv.inp file. See mk_tPDB.tcl"}
print "CFACV) nCV $nCV"
if {[info exists TAMDverbose]} {print "CFACV) cvList: $cvList"}


# Allocate dataspace in C code
set ds [NewDataSpace $nCntr $nCV $restr(num) $seed]

# Add restraints to the dataspace structure and setup restraint global variable
rlist_setup $nCV

# Add CVs to the dataspace structure
cvs_setup

# Add chapeau to the dataspace structure and setup chapeau global variable
# TODO: improve chapeau global variable
if {[info exists CFACV_doAnalyticalCalc]} {
  if {$CFACV_doAnalyticalCalc == 1} {
    if {![info exists USETAMDFORCES]} {set USETAMDFORCES 0}

    # Replica exechange mode: 
    if {![info exists NUMREP]} {set NUMREP 1}
    # TODO: NUMREP podria ser un numero de chapeaus distintos (estilo el
    # problema sodio cloro) a obtener en una simulacion. En ese caso
    # corregir.

    # Saving for allocate chapeau functions
    chapeau_setup $NUMREP $ds $XSCFILE $SPLINEMIN $NKNOTS $CUTOFF $BEGINEVOLVEPARAMETERS $USETAMDFORCES $BINREPORTPARAMFILE $BINREPORTPARAMFREQ $BINOUTPUTLEVEL $LAMUPDATEINTERVAL
  }
}

# read z values from a restart file
set first 1
if {[info exists restart_root]} {
  cfacv_loadstate $ds $NUMREP $restart_root
  set first 0
}

# Some default values regarding output
if {![info exists TAMDof]}          {set TAMDof 1 }
if {![info exists TAMDbinof]}       {set TAMDbinof 1 }
if {![info exists TAMDoutputlevel]} {set TAMDoutputlevel 3}

set TAMDoutputFileFP 0
if {[info exists TAMDoutputFile]} {
    set TAMDoutputFileFP [my_fopen $TAMDoutputFile "w"]
    print "CFACV) Opened TAMD output file $TAMDoutputFile"
} else {
    set TAMDoutputFileFP [my_fopen stdout "w"]
    print "CFACV) TAMD output to stdout"
}

set TAMDbinOutputFileFP 0
if {[info exists TAMDbinOutputFile]} {
    set TAMDbinOutputFileFP [my_binfopen $TAMDbinOutputFile "w" $TAMDoutputlevel $ds]
    print"CFACV) Binary TAMD output to $TAMDbinOutputFile"
}

# define the "calcforces" function which is called 
# at each MD timestep
proc calcforces { } {

    global ds
    global groups
    global first
    global TAMDof
    global TAMDbinof
    global TAMDoutputlevel
    global TAMDoutputFileFP
    global TAMDbinOutputFile
    global TAMDbinOutputFileFP

    # load coordinates of requested atoms into associative array
    # this is a NAMD builtin
    loadcoords p
    # print $p(g1)

    # perform the update that transmits forces
    Tcl_UpdateDataSpace $ds p $groups $first [getstep]
    if {$first==1} { set first 0 }

    # report if requested
    if {[expr {[getstep]%$TAMDof == 0}]} {
      DataSpace_ReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDoutputFileFP
    }

    # report if requested
    if {[info exists TAMDbinOutputFile] && [expr {[getstep]%$TAMDbinof == 0}]} {
      DataSpace_BinaryReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDbinOutputFileFP
    }
}

