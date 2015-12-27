
proc chapeau_setup { NUMREP ds XSCFILE min N max begin_evolve begin_solve usetamdforces outputfile outfreq outlevel nupdate} {
    global chape

    # FIXME: This code was here to allow compute chapeau functions separatedly
    # for different pair types of particles. For instance, this allow to
    # recover SOD SOD, CLA CLA and SOD CLA pair potentials in 1 TAMD
    # simulation. Each index has a number in ch_id which allow to sort the pair
    # in the different chapeau objects on the c code.  From the studies with
    # SOD CLA, this pair potentials will be OK only if the ficticius
    # temperature is the same that the real one.  On the other hand, a better
    # way to achive this is needed (without saving a lot of numbers in ch_id).
    # For understand how this worked, see the previous versions of the code.
    #
    # # chplist holds the list of pair potential types 
    # # e.g. {{SOD CLA} {CLA CLA} {SOD SOD}}
    # if {![info exists chlist]} {set chlist "{}"}
    # set chnum [llength $chlist]
    #
    # # Pair potentinal interaction
    # intListToArray_Data [DataSpace_chid $ds ] $ch_id
    # print "CFACV) The pair chapeau (map in vector) is:"
    # for {set i 0} { $i < $nCntr } { incr i } {
    #   set aux ""
    #   for {set j 0} { $j <= $i } { incr j } {set aux "$aux  "}
    #   for {set j [expr $i+1]} { $j < $nCntr } { incr j } {
    #     set aux "$aux [lindex $ch_id [expr $j+($nCntr-2)*$i-($i-1)*$i/2-1]]"
    #   }
    #   puts $aux
    # }

    # FIXME: pairwise analytical potential
    # DataSpace_SetupPairCalc $ds $min $N $coff $nlcoff $begin_evolve $usetamdforces $outputfile $outfreq $outlevel $nupdate $chnum

    # TODO: I guess one setup subroutine per chapeau object would be nice, so each chapeau can have different grid
    # and be asociated with differen CV sets

    # Currently only option is a 1D grid
    DataSpace_Setup1Dchapeau $ds $NUMREP $min $N $max $begin_evolve $begin_solve $usetamdforces $outputfile $outfreq $outlevel $nupdate

    # Get the adress of the chapeaus
    for {set i 0} {$i<$NUMREP} {incr i} {
      set chape($i.address) [DataSpace_get_chapeauadress $ds $i]
    }

    if [string equal $XSCFILE "off"] {
      print "CFACV) DEBUG: Tcl_InitializePairCalc nobox"
      DataSpace_SetupPBC $ds 0  0 0 0  0 0 0
    } else {
      set LL [my_getcellsize $XSCFILE]
      set O [my_getorigin $XSCFILE]
      print "CFACV) DEBUG: Tcl_InitializePairCalc box size [lindex $LL 0] [lindex $LL 1] [lindex $LL 2]"
      print "CFACV) DEBUG: Tcl_InitializePairCalc origin   [lindex $O 0] [lindex $O 1] [lindex $O 2]"
      DataSpace_SetupPBC $ds 1 [lindex $O 0] [lindex $O 1] [lindex $O 2] [lindex $LL 0] [lindex $LL 1] [lindex $LL 2] 
    }
}

