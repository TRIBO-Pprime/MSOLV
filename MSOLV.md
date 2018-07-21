---
project: MSOLV-fortran
summary: MSOLV-fortran -- A modern Fortran API for Sparse Direct Solvers, as part of <br/><br/> ![MUSST_img](media/MUSST_long_petit.png)
author:  Arthur Francisco - Noël Brunetière
email:   arthur.francisco@univ-poitiers.fr
website: https://www.pprime.fr/?q=fr/recherche-scientifique/d3/mecanique-des-interfaces-lubrifiees
github: https://github.com/Arthur-Francisco
project_github: https://github.com/TRIBO-Pprime/MSOLV
include:    ./inc
src_dir:    ./src
exclude: hsl_common90.f90
         hsl_ddeps90.f90
         hsl_ma48d.f90
output_dir: ./docs
media_dir:  ./img
favicon:    ./img/logo.png
docmark:        !
docmark_alt:    #
predocmark:     >
predocmark_alt: <
display: public
         protected
         private
source: true
graph:  true
search: true
sort: src
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc(anchorlink=False)
               markdown.extensions.smarty(smart_quotes=False)
---

-----------------

[TOC]

Brief description
-----------------

A generic API for using Sparse Direct Solvers, written in modern Fortran. It is designed to work with:

* [MUMPS](http://mumps.enseeiht.fr/index.php?page=doc) - **MU**ltifrontal **M**assively **P**arallel sparse direct **S**olver [^1]
* [UMFPACK](http://faculty.cse.tamu.edu/davis/suitesparse.html) - **U**nsymmetric **M**ulti**F**rontal method solver [^2]
* [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu) - LU decomposition with partial pivoting and triangular system solves through forward and back substitution [^3]
* [MA48](http://www.hsl.rl.ac.uk/catalogue/ma48.html) - Sparse unsymmetric system solver [^4]

Requested libraries
-------------------

* *Generic*</br>
	+ ```librefblas.a```
	+ ```liblapack.a```</br></br>
* **MUMPS**</br>
	+ ```libmetis.a```
	+ ```libmpiseq.a```
	+ ```libmumps_common.a```
	+ ```libdmumps.a```</br></br>
* **UMFPACK**</br>
	+ ```libamd.a```
	+ ```libcamd.a```
	+ ```libcolamd.a```
	+ ```libccolamd.a```
	+ ```libcholmod.a```
	+ ```libsuitesparseconfig.a```
	+ ```libumfpack.a```</br></br>
* **SuperLU**</br>
	+ ```libsuperlu.a```</br></br>
* **MA48**</br>
	+ ```libhsl_ma48.a```

Example of use
--------------

[[test_solvers]] proposes to test the solvers on two systems, a very small one and a "bigger" one. Some solvers have an "analyze" step during which symbolic factorization is performed.
Hence if the sparsity of a new system does not change much, the symbolic factorisation is reused for the numerical factorization. To assess the time gained, after the first resolution other
resolution are proposed with the same matrix pattern.
After each resolution, the error is assessed remultipying the solution by the matrix and comparing it to the right-hand side.

Data format
-----------

By default, 

* **UMFPACK** and **SuperLU** are configured to read **C**ompressed **C**olumn [Harwell/Boeing](http://netlib.org/linalg/html_templates/node92.html) formatted data,</br></br>
* **MUMPS** reads elemental CC formatted data,</br></br>
* **MA48** reads sparse triplet formatted data.

Limitations
-----------

The package **MUSST** is developped to handle some lubrication problems thanks to a multiscale approach. During the processe, sparse linear systems are solved either on concurrent threads or using all threads. In the first case, 'OpenMP' is not usefull because each thread is dedicated to the resolution of its own bunch of systems (**B**ottom **S**cale systems). In the second case, 'OpenMP' is desirable to speed up the resolution (**T**op **S**cale systems). For the moment, **MUSST** isn't designed for distributed shared memory architectures. 
Therefore the following choices are made:

* **MUMPS** is built with its *sequential* configuration (*OpenMP*, but not *MPI*). It is dedicated to **TS** systems.</br></br>
* **BLAS** is built without *OpenMP*</br></br>
* **UMFPACK**, **SuperLU** and **MA48** are built without *OpenMP*. They are dedicated to **BS** systems.

Prerequesites
-------------
The solver libraries must be downloaded and statically compiled together with their prerequisites:

* **UMFPACK**</br>
	*BLAS*, *LAPACK*</br></br>
* **MUMPS**</br>
	*BLAS*, *SCOTCH*</br></br>
* **SuperLU**</br>
	*BLAS*</br></br>
* **MA48**</br>
	*BLAS*
	
Wrappers
--------

* **UMFPACK** </br>
	A module [[mumfpack]] has been developped by [Ladislav Hanyk](http://geo.mff.cuni.cz/~lh/Fortran/UMFPACK/). It has been checked with UMFPACK 5.6.2, but it also works fine here with UMFPACK 5.7.6.</br></br>
* **MUMPS**</br>
	A module [[mumps_wrapper]] is created with the content of ```dmumps_struc.h``` provided by MUMPS and it includes ```mpif.h```, also provided by MUMPS. It is to be noticed that almost all global variables are deactivated in ```mpif.h```.</br></br>
* **SuperLU**</br>
	A module [[sulu_wrapper]] is created from scratch because following SuperLU fortran examples (with the wrapper provided by SuperLU), we were not able to continuously solve systems with systematic memory release. As a consequence, the memory needed by MUSST increased untill the programs would stop.</br></br>
* **MA48**</br>
	In its *double* flavour, *HSL_MA48* needs 3 files: ```common90.f90```, ```ddeps90.f90``` and ```hsl_ma48d.f90```, renamed [[hsl_common90.f90]] and [[hsl_ddeps90.f90]] for the two first.
	In order to work with arrays independently from the solvers, the global solver type [[MAT_SOLV]] contains the system arrays --like ```eltptr```, ```eltvar```, ```a_elt```, etc.-- pointed by a solver type, for instance [[ZD11_TYPE]]. Therefore, in version 3.3 of *HSL_MA48*, the attribute ```allocatable``` (in ZD11_TYPE) of ```row```, ```col```, ```ptr``` and ```val``` is changed to ```pointer```.

Make
----
General syntax:
```bash
./comp.sh ALL:s/a FORT:gfortran-x DEBUG:yes/no GPROF:yes/no
```
Examples with gfortran-7:

* first rebuild all, then make the 'debug' executable ```prg``` for profiling
```bash
./comp.sh a gfortran-7 yes yes
```

* just build the modified sources, then make a productive executable
```bash
./comp.sh s gfortran-7 no no
```

Full description
-----------------

The aim of the present project is to provide generic subroutines for the direct full resolution of sparse linear systems:

* solver initialization
```fortran
call solve_syst(mat=system_type, step='ini')
```
* system analyze
```fortran
call solve_syst(mat=system_type, step='ana')
```
* system factorization
```fortran
call solve_syst(mat=system_type, step='fac')
```
* system solution
```fortran
call solve_syst(mat=system_type, step='sol')
```
* memory release
```fortran
call solve_syst(mat=system_type, step='fre')
```
* solver finalization
```fortran
call solve_syst(mat=system_type, step='end')
```

```system_type``` is of complex type [[MAT_SOLV]] built to deal with *MUMPS*, *UMFPACK*, *SuperLU* and *MA48* requirements.


License
-------



[^1]:
	CeCILL-C license
[^2]:
	GNU GPL
[^3]:
	New BSD-3
[^4]:
	HSL software is strictly intended for Personal academic use on the download page


