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

set BIAS 0
enabletotalforces

proc calcforces { } {
  # calcforces is call each timestep before move the system, and 1
  # aditional time after the system has moved.

  global ds
  global groups
  global first
  global TAMDof
  global TAMDbinof
  global TAMDoutputlevel
  global TAMDoutputFileFP
  global TAMDbinOutputFile
  global TAMDbinOutputFileFP
  global BIAS
         
  # Segun mis test y deducciones en el primer paso namd hace:
  # 1- llama a calcforces
  # 2- calcula las energías (hasta aca un un `run 0`)
  # 3- mueve
  # 4- llama a calcforces
  # 5- vuelve a clacular las energías e imprimendolas en pantalla.
  # 6- ejecuta el callback
  #
  # Hasta aca coincide la posicion con la energia
  # en pantalla, y hasta aca 2 llamados a
  # calcforces y luego el callback.
  #
  # El segundo paso evita repetir el paso anterior, y luego:
  # 1- llama a calcforces
  # 2- mueve
  # 3- vuelve a calcular las energías
  # 4- ejecuta el callback
  #
  # Notar que recien aca en 1 coincide posicion con la energía guardad
  # en el callback del paso anterior. Por eso conviene si se usa el
  # calllback (como en una dinamica acelerada) evitar los primeros 2 pasos
  # de registrar nada.

  # load coordinates of requested atoms into associative array
  # this is a NAMD builtin
  loadcoords p
  loadforces fur
  loadtotalforces f
  # puts "[getstep]"
  # if {$first==0} {
  #   # puts "FUERZAS3 $f"
  #   # loadtotalforces f
  #   puts "FUERZAS3 $f(g1)"
  #   # puts "FUERZAS3 $fur(g1)"
  # #   if {[getstep]>2} exit
  # }


  # print $p(g1)

  # perform the update that transmits forces

  Tcl_UpdateDataSpace $ds p $groups $first [getstep] $BIAS
  if {$first==1} { set first 0 }

  # report if requested
  if {[info exists TAMDbinOutputFile] && [expr {[getstep]%$TAMDbinof == 0}]} {
    DataSpace_BinaryReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDbinOutputFileFP
  }
}

