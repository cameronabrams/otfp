# Scripts to be call by NAMD at each MD timestep.
# The initialization script is a separated file in order to execute any desired
# command after initialize and before run.

# Set up the centers as "groups". This give the list groups with the groups ids
# and only can be done in the script call by tclforces.
set groups {}
for {set i 1} { $i <=$pdb(nctrs) } { incr i } {

  # addgroup is a NAMD command that return a "gN" id with N a small integer
  set gid [addgroup $pdb($i.atoms)]
  lappend groups $gid

  if [info exists TAMDverbose] {print "addgroup $i  $pdb($i.atoms)"}
}
 

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
  if {[info exists TAMDbinOutputFile] && [expr {[getstep]%$TAMDbinof == 0}]} {
    DataSpace_BinaryReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDbinOutputFileFP
  }
}
 
