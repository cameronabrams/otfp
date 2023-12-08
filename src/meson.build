cc = meson.get_compiler('c')
m = cc.find_library('m', required: false)

lapack = dependency('lapack')
blas = dependency('blas')
tcl = dependency('tcl')

swig = find_program('swig')
swig_gen = generator(
  swig,
  output: ['@BASENAME@_wrap.c'],
  arguments: ['-tcl8', '-o', './@OUTPUT@','@INPUT@'],
  depfile: '@BASENAME@.d',
)

install_headers('calcforces.tcl')
install_headers('cfacv.tcl')
install_headers('chapeau.tcl')
install_headers('cvs.tcl')
install_headers('data.tcl')
install_headers('restraints.tcl')

 

cfacv_src = files(
  'measurements.c',
  'measurements.h',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  'centers.c',
  'centers.h',
  'cfacv.c',
  'cfacv.h',
  'cfacv_obj.h',
  'cvs.c',
  'cvs.h',
  'cvs_obj.h',
  )

cfacv_src += swig_gen.process('cfacv.i')

shared_library(
  'cfacv',
  cfacv_src,
  dependencies: [lapack, blas, m, tcl],
  c_args : '-DVERSION="@0@"'.format(meson.project_version()),
  install: true,
  install_dir: 'lib/',
)



cvs_src = files(
  'measurements.c',
  'measurements.h',
  'cvs.c',
  'cvs.h',
  'cvs_obj.h',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
)

cvs_src += swig_gen.process('cvs.i')

shared_library(
  'cvs',
  cvs_src,
  dependencies: [lapack, blas, m, tcl],
  c_args : '-DVERSION="@0@"'.format(meson.project_version()),
  install: true,
  install_dir: 'lib/',
)  



catbinsp_src = files('catbinsp.c')
 
executable('catbinsp', 
  catbinsp_src, 
  dependencies: [m],
  install: true,
)
 


chapsim_src = files(
  'simetrize.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
)


executable('chapsim', 
  chapsim_src, 
  dependencies: [lapack, blas, m],
  install: true,
)
 


chapconv_src = files(
  'convchapeau.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  'chapeau_old.h',
  ) 

executable('chapconv', 
  chapconv_src, 
  dependencies: [lapack, blas, m],
  install: true,
  )



chapadd_src = files(
  'addchapeaus.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  )

executable('chapadd', 
  chapadd_src, 
  dependencies: [lapack, blas, m],
  install: true,
  )



chapcrop_src = files(
  'cropchapeaus.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  )

executable('chapcrop', 
  chapcrop_src, 
  dependencies: [lapack, blas, m],
  install: true,
  )



