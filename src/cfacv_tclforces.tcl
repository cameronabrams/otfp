# tclforces file for CFACV
# Cameron F Abrams
# 2009-14

#### Check for necessary parameters set in main config. file

# # Runing alone
# set CFACV_BASEDIR $env(PWD)
# set DIRS_INPUT ..
# set labelPDB $DIRS_INPUT/label.pdb 
# set cvINP $DIRS_INPUT/cv.inp
# proc print { arg } {
#   puts $arg
# }
# set TAMDof                  10
# set TAMDoutputlevel         0
# set CFACV_doAnalyticalCalc  1
# set XSCFILE                 $DIRS_INPUT/tip3pE.xsc
# set CUTOFF                  7.0
# set NLCUTOFF                8.0
# set BEGINEVOLVEPARAMETERS   99
# set REPORTPARAMFREQ         100
# set SPLINEMIN               0.0
# set NKNOTS                  141
# set BINREPORTPARAMFREQ      100
# set BINREPORTPARAMFILE      $DIRS_INPUT/alone.bsp
# set BINOUTPUTLEVEL          3
# set initKnotsINP            $DIRS_INPUT/wcaknots
 
 
# # Runing from NAMD
# set CFACV_BASEDIR $PWD/OTFP
 
# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/cfacv.tcl

cfacv_banner NAMD

set tripped 0
foreach key {labelPDB cvINP} {
    if {![info exists $key]} {
	print "CFACV) ERROR: you must set $key in the NAMD config. file."
	exit
    }
}

# pairmask holds the list of pair potential types 
# e.g. {{SOD CLA} {CLA CLA} {SOD SOD}}
if {![info exists pairmask]} {set pairmask ""}


#### Get the groups using the addgroup of tclforces

set serArray {}; # must have for addgroup
set masses {}
set ch_id {}

# Read the template PDB file that identifies subdomain memberships
set nCntr [read_centersPDB $labelPDB serArray masses pairmask ch_id]
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

# Set up list of CV's
set cvList {}
set nCV [read_cvs $cvINP cvList $nCntr]
print "CFACV) nCV $nCV"

#Now, cvList is some like {CARTESIAN_X 0} {CARTESIAN_Y 0} ....
if {[info exists TAMDverbose]} {print "CFACV) cvList: $cvList"}

if {!$nCV} {
    print "CFACV) ERROR: Perhaps you need to use mk_tPDB.tcl to generate the cv.inp file?"
}

# Set up list of restraints
set rList {}
if {[info exists restrINP]} {
    set nR [read_restraints $restrINP $nCV rList]
    print "CFACV) $restrINP : nRestraints $nR"
} else {
    set nR [create_single_cv_restraints $nCV rList $restrPARAMS]
    print "CFACV) $nR single-cv restraints created:"
}

if {[info exists TAMDverbose]} {foreach r $rList {print "CFACV) $r"}}
 
#Now, rListis something like "{{1 0 0 0 0} {k 1600} {AMDkT 0.19} {TAMDgamma
#1} {TAMDdt 0.002}} {{0 1 0 0 0} {k 1600} {AMDkT 0.19} {TAMDgamma 1}
#{TAMDdt 0.002}} ... "


# set the pseudorandom number generator seed
#if {![info exists seed]} {set seed [clock clicks]}
print "CFACV) setting seed to $seed"
   
# declare and allocate data space
# Here all the previous information is stacked
set ds [Tcl_NewDataSpace $nCntr $cvList $rList $seed]

# if intercenter pair calcs are needed
if {[info exists CFACV_doAnalyticalCalc] && $CFACV_doAnalyticalCalc == 1} {

    if {![info exists USETAMDFORCES]} {set USETAMDFORCES 0}

    # currently only option is a pairwise analytical potential
    set aux  [llength $pairmask]
    Tcl_InitializePairCalc $ds $XSCFILE $CUTOFF $NLCUTOFF $BEGINEVOLVEPARAMETERS $USETAMDFORCES $REPORTPARAMFREQ $SPLINEMIN $NKNOTS $BINREPORTPARAMFILE $BINREPORTPARAMFREQ $BINOUTPUTLEVEL $LAMUPDATEINTERVAL $aux

    # FIXME. Add a index to load a  initial knots file for each chapeau
    # if {[info exists initKnotsINP]} {Tcl_DataSpace_InitKnots $ds $initKnotsINP}

    # Pair potentinal interaction
    intListToArray_Data [DataSpace_chid $ds ] $ch_id
    print "CFACV) The pair chapeau (map in vector) is:"
    for {set i 0} { $i < $nCntr } { incr i } {
      set aux ""
      for {set j 0} { $j <= $i } { incr j } {set aux "$aux  "}
      for {set j [expr $i+1]} { $j < $nCntr } { incr j } {
        set aux "$aux [lindex $ch_id [expr $j+($nCntr-2)*$i-($i-1)*$i/2-1]]"
      }
      puts $aux
    }
              
}

# read z values from a restart file
set first 1
if {[info exists restartINP]} {set first [Tcl_Reinitialize $ds $restartINP]}

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
