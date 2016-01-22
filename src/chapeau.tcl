
proc chapeau_setup { NUMREP ds DIMEN min N max begin_evolve begin_solve usetamdforces outputfile outfreq outlevel nupdate} {
    global chape

    # Arrays
    set a1 [ListToArray $min] 
    set a2 [intListToArray $N] 
    set a3 [ListToArray $max] 

    # Currently only option is a 1D grid
    DataSpace_SetupChapeau $ds $NUMREP $DIMEN $a1 $a2 $a3 $begin_evolve $begin_solve $usetamdforces $outputfile $outfreq $outlevel $nupdate

    delete_array $a1
    delete_arrayint $a2
    delete_array $a3

    # Get the adress of the chapeaus
    for {set i 0} {$i<$NUMREP} {incr i} {
      set chape($i.address) [DataSpace_get_chapeauadress $ds $i]
    }

}

