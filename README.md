This is a [fork](https://github.com/alexispaz/otfp) from the original OTFP code
hosted in Cameron Abrams' [github](https://github.com/cameronabrams/otfp).

# On the fly parameterization 

OTFP computes free-energy profiles in MD simulations via temperature-acceleration. 
 
It is written in C and use [SWIG](http://www.swig.org) to provida a TCL
interface to be used with [NAMD](www.ks.uiuc.edu/Research/namd) via TCLforce.
 
## Build system

There is a `src/makefile` file that can be used to build libraries and
executables. It is also posible to compile using meson, e.g.:

    git clone https://github.com/alexispaz/otfp.git
    cd otfp
    meson setup build 
    cd build 
    meson compile -v 
    meson install 
    cd ..

Default prefix is `/usr` but can be changed with `--prefix`, for instance:

    meson setup --prefix=$(realpath ./usr) build 

To run the examples, set the `CFACV_BASEDIR` environment variable to the OTFP
prefix:

    cd examples/OTFP/
    export CFACV_BASEDIR=/usr
    ./test.sh

## Citations

The OTFP method is explained in detail in the original publication:

- Abrams and Vanden-Eijnden. 
_"On-the-fly free energy parameterization via temperature accelerated molecular
dynamics."_
Chem. Phys. Lett., 547 (2012) 114–119. https://doi.org/10.1016/j.cplett.2012.07.064

The extension of OTFP in 2D is explained in:

- Paz, Maragliano, and Abrams. 
_"Effect of Intercalated Water on Potassium Ion Transport through Kv1.2
Channels Studied via On-the-Fly Free-Energy Parametrization."_
J. Chem. Theory Comput., 14(5) (2018) 2743–2750. https://doi.org/10.1021/acs.jctc.8b00024
 
The combination of OTFP with Replica Exchanges is described in:

- Paz, Vanden-Eijnden and Abrams. 
_"Polymorphism at 129 dictates metastable conformations of the human prion
protein N-terminal β-sheet."_ 
Chem. Science, 8(2) (2016) 1225–1232. https://doi.org/10.1039/c6sc03275c
 
Other publications with OTFP method are:

- Alberini, Paz, Corradi, Abrams, Benfenati and Maragliano.
"Molecular Dynamics Simulations of Ion Permeation in human NaV Channels."
J. Chem. Theory Comput., 19(10) (2023) 2953–2972. https://doi.org/10.1021/acs.jctc.2c00990

- Paz and Abrams. 
_"Free Energy and Hidden Barriers of the β-Sheet Structure of Prion Protein."_
J. Chem. Theory Comput., 11(10) (2015) 5024–5034. https://doi.org/10.1021/acs.jctc.5b00576

- Paz and Abrams. 
_"Testing Convergence of Different Free-Energy Methods in a Simple Analytical
System with Hidden Barriers."_
Computation, 6(2) (2018) 27.  https://doi.org/10.3390/computation6020027

## License

It mantains the same BSD-2-Clause license of the original code at [repo](https://github.com/cameronabrams/otfp).

Original code was developed with support of the National Science Foundation
through grant DMR-1207389.
 
## Known issues

In [SWIG documentation](http://www.swig.org/Doc1.3/Tcl.html) it is written: 

> "As a final complication, a major weakness of C++ is that it does not define any
> sort of standard for binary linking of libraries. This means that C++ code
> compiled by different compilers will not link together properly as libraries
> nor is the memory layout of classes and data structures implemented in any kind
> of portable manner. In a monolithic C++ program, this problem may be unnoticed.
> However, in Tcl, it is possible for different extension modules to be compiled
> with different C++ compilers. As long as these modules are self-contained, this
> probably won't matter. However, if these modules start sharing data, you will
> need to take steps to avoid segmentation faults and other erratic program
> behavior. If working with lots of software components, you might want to
> investigate using a more formal standard such as COM."

Thus, it is recommended to compile OTFP using the same compiler used to compile
NAMD. 

Another problem that can arise is a segmentation fault when the TCL script of
NAMD input loads the compiled `cfacv.so` shared library. This might be related with
differences in the linking of the tcl library in NAMD and in the `cfacv.so`
library. To find which tcl library is using NAMD execute `ldd ./namd2`.  Use
that library to compile `cfacv.so`.

If the tcl library is not found by ldd, it is probably linked static. In that
case `nm ./namd2 | grep Tcl` will show the common symbols defined in the tcl
library (e.g. `Tcl_AppendResult`). This make difficult to known which tcl
version was used to compile NAMD. If using the tcl library of the actual system
gives segmentation fault, sometimes the issue disappears if the user compile
its own TCL library from source and use that one to compile of `cfacv.so`.

## Examples of additional tools to handle chapeau object files

Convert old format of chapeau object files, create `converted.ch`

    ./chapconv run12.ch0 

Cropping or expanding the CV space, create `croped.ch`

    ./chapcrop cryo.ch 4 18.40076 4 18.40076

Attempt to add chapeau object files

    ./chapadd croped_*ch

If attempt fail, it might be needed to adjust the size of chapeau objects to
have the same number of knots using `chapcrop`




