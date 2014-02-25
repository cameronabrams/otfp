# tclforces file for CFACV
# Cameron F Abrams
# 2009-14

# check for base directory name environment variable;
# if not set, use default
#if {[info exists env(CFACV_BASEDIR)]} {
#    set CFACV_BASEDIR $env(CFACV_BASEDIR)
#} else {
    set HOME $env(HOME)
    set CFACV_BASEDIR ${HOME}/cfacv_otfp
#}

# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/cfacv.tcl

cfacv_banner NAMD

# Check for necessary parameters set in main config. file
set tripped 0
foreach key {labelPDB cvINP} {
    if {![info exists $key]} {
	print "CFACV) ERROR: you must set $key in the NAMD config. file."
	set tripped 1
    }
}
if {$tripped} {
    exit
}

set serArray {}; # must have for addgroup
set masses {}

# read the template PDB file that identifies subdomain memberships
set nCntr [read_centersPDB $labelPDB serArray masses]
print "CFACV) nCenters $nCntr  masses $masses"

# Set up the subdomains as "groups" for tclforces
set groups {}
for {set i 0} { $i < $nCntr } { incr i } {
    if {[info exists TAMDverbose]} {
	print "addgroup $i :"
	print "   [lindex $serArray $i]"
    }
    set gid [addgroup [lindex $serArray $i]]
    lappend groups $gid
}

# Set up list of CV's
set cvList {}
set nCV [read_cvs $cvINP cvList $nCntr]
print "CFACV) nCV $nCV"
if {[info exists TAMDverbose]} {
    print "CFACV) cvList: $cvList"
}

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

if {[info exists TAMDverbose]} {
    foreach r $rList {
	print "CFACV) $r"
    }
}

# set the pseudorandom number generator seed
#if {![info exists seed]} {
#    set seed [clock clicks]   
    print "CFACV) setting seed to $seed"
#}

# declare and allocate data space
set ds [Tcl_NewDataSpace $nCntr $cvList $rList $seed]

# if intercenter pair calcs are needed
if {[info exists CFACV_doAnalyticalCalc] && $CFACV_doAnalyticalCalc == 1} {
    if {![info exists USETAMDFORCES]} {
	set USETAMDFORCES 0
    }
    # currently only option is a pairwise analytical potential
    Tcl_InitializePairCalc $ds $XSCFILE $CUTOFF $NLCUTOFF $BEGINEVOLVEPARAMETERS $USETAMDFORCES $REPORTPARAMFREQ $SPLINEMIN $NKNOTS $BINREPORTPARAMFILE $BINREPORTPARAMFREQ $BINOUTPUTLEVEL $LAMUPDATEINTERVAL
    if {[info exists initKnotsINP]} {
	Tcl_DataSpace_InitKnots $ds $initKnotsINP
    }
}

# read z values from a restart file
set first 1
if {[info exists restartINP]} {
    set first [Tcl_Reinitialize $ds $restartINP]
}

if {[info exists TAMDof]} {
    set reportFreq $TAMDof
} else {
    set reportFreq 1
}

if {[info exists TAMDbinof]} {
    set binReportFreq $TAMDbinof
} else {
    set binReportFreq 1
}


if {![info exists TAMDoutputlevel]} {
    set TAMDoutputlevel 3; # default output Z and Theta
}

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
    global reportFreq
    global binReportFreq
    global TAMDoutputlevel
    global TAMDoutputFileFP
    global TAMDbinOutputFile
    global TAMDbinOutputFileFP

    # load coordinates of requested atoms into associative array
    loadcoords p

    # perform the update that transmits forces
    Tcl_UpdateDataSpace $ds p $groups $first [getstep]
    if {$first==1} { set first 0 }

    # report if requested
    if {[expr {[getstep]%$reportFreq == 0}]} {
	DataSpace_ReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDoutputFileFP
    }

    # report if requested
    if {[info exists TAMDbinOutputFile] && [expr {[getstep]%$binReportFreq == 0}]} {
	DataSpace_BinaryReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDbinOutputFileFP
    }
}
