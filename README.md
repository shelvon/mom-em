[GNU GENERAL PUBLIC LICENSE Version 2, June 1991](./LICENSE.md)

# USER'S GUIDE

The project name is changed to **[mom-em](https://github.com/shelvon/mom-em)** (which stands for **the Method of Moments for Electromagnetics**)

This project is intended to extend functionalities of the previous **[for90-mom2](https://github.com/jmakitalo/for90-mom2)** project that was created by previous group member Jouni Mäkitalo. The current author Xiaorun ZANG (shelvon) joined, as a postdoctoral researcher, in the Nonlinear Optics Group at Tampere University of Technology in 2016, since when he is in charge of the management and extension of this open source code.

**Current author information**

Xiaorun (shelvon) ZANG (xiaorun.zang@tuni.fi)  
Postdoctoral Researcher at Tampere University (2016-)  
Photonics Laboratory, Physics Unit, Tampere University,  
P.O. Box 692, FI-33014 Tampere, Finland  

**Original author information**

Jouni Mäkitalo (jmakitalo15@gmail.com)  
Researcher at Tampere University of Technology (2011-2015)  
Department of Physics  
Optics Laboratory

**References**

The code was initially developed in jouni's Ph.D. thesis:

"Boundary Integral Operators in Linear and Second-order Nonlinear Nano-optics"  
[http://URN.fi/URN:ISBN:978-952-15-3539-0](http://URN.fi/URN:ISBN:978-952-15-3539-0)

**Disclaimer**

*Although the program has been used to successfully model several cases and the results compare well with the results for a spherical particle, the authors do not take any responsibility for the correctness of results obtained by the end-user.*

**Preface**

This is the user's guide to the Method of Moments in Electromagnetics solver implemented in [Fortran](https://www.wikiwand.com/en/Fortran).

This program is based on the theory and expressions derived in *Master of Science Thesis (Jouni Mäkitalo, 2011)*.  
The nonlinear surface response code was developed in his *Ph.D. Thesis (Jouni Mäkitalo, 2015)*.

The program is made for solving linear and weakly nonlinear scattering problems, where the particle size is on the order of the wavelength and whose material is described by a piece-wise constant complex permittivity. The following list gives the main features and application targets:

- Computation of full-wave solutions to electromagnetic scattering problems
- PMCHWT boundary element formulation of dielectric and lossy piece-wise homogeneous scatterers
- Deduce scattering and absorption cross-sections, near-fields and radiation patterns
- Excitations: plane-wave, focused Gaussian beams, electric dipole, general vector beams of paraxial/nonparaxial type
- Enforce symmetries and decompose solutions from their irreducible representations
- Second-order nonlinear scattering with surface electric dipolar source
- Second/third-order nonlinear scattering with bulk dipolar source
- Multidomain problems: a domain for each constant permittivity, the domains may
  be nested and may be touching
- Problems with 2-D periodicity in the xy-plane, deduce near-fields, reflectivity and transmittance
- Solve a problem for multiple sources quickly by reusing matrix factorization.
- Nonlinear beamscanning: deduce the reflected second-harmonic power when a focused beam is
  scanned over a rectangular domain
- Eigenmodes of flat PEC scatterers/antennas.
- Volume integral formulation of scattering from inhomogeneous objects. Uses SWG basis and Duffy transform.

## Building

*Fortran compilers* (the following two compilers under Linux/macOS environment are tested):  
- `ifort`  
The code has been successfully built with the Intel Fortran compiler of different versions (ifort 12.1.2-18.0.5).  
- `gfortran`  
The compatibility with gfortran has been improved later, currently the code should be complied under gfortran 4.8.5 (or later).

If the end-user runs into errors, he/she is welcomed to contact us.

*Libraries* that should be installed separately:  
- `HDF5`, the fortran support must be enabled. One also needs to set `HDF5_ROOT` environment variable to the `HDF5` installation folder.  
- `LAPACK`, during the development, *Intel's Math Kernel Library (MKL)* packaging of `LAPACK` was used. One needs to set `MKLROOT` environment variable to the `MKL` installation folder.  
- `OPENMP`, support for parallel computation. This library is generally  included in compilers.  

*Libraries* that have been included in this project:  
- `fson`, the *Fortran 95 JSON Parser* is placed under `src/fson` forked from [https://github.com/josephalevin/fson](https://github.com/josephalevin/fson).  
- `toms644`, computes all major **Bessel functions of a complex argument and nonnegative order**.  ALGORITHM 644, COLLECTED ALGORITHMS FROM ACM. THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 21, NO. 4, December, 1995, P.  388--393. [http://www.netlib.org/toms-2014-06-10/644](http://www.netlib.org/toms-2014-06-10/644)  

*Build steps*:

```
mkdir bin
mkdir build
cd build
FC=ifort cmake ../src/
make
```
Alternatively, set `FC=gfortran` to use gfortran compiler, and use `make -j -k` to speed up the compilation.

The executeable named `mom` will be placed in `bin/` directory.

## Basics usage

### Parallel computation

The program supports symmetric multiprocessing via OpenMP (when installed). To get n threads, set the environment variable as follows.

	"OMP_NUM_THREADS=n"
### Read parameters

There are two ways to use this program.
1. Read parameters from a json file.  
	The program will treat the json file as a structure type variable, and read all the necessary nodes, subnodes, ..., recursively.    
		`./mom input.json`  

2. Get commands line by line  
		- either interactively  
			`./mom`  
		The program will display a prompt and will receive commands.  
		- or redirection from an input file.  
			`./mom < input`  
		In Unix-like operating systems you can pass the commands via a text file with stdin redirection.

#### 1. Read parameters from JSON file (Recommended)

There are several advantages of using this new method for reading parameters.  
- The order of nodes in a json file doesn't matter.  
- It is easy to use other glue languages, such as bash shell (requires jq, [https://stedolan.github.io/jq/](https://stedolan.github.io/jq/)), to control the generation or modification of a json file thus to sweep some parameters.  
- The data in a json file has a tree structure and it is more readable.  

A typical json root object looks like as follow.  
(Each *"node"* is preceded by a *"#node"* which describes the *"node"*.) 

```javascript
{
  "#name": "[required] string type (maximum length, 120 characters) that specifies the job's name",
  "name": "jobname",
  "#geom": "[required if .simulation.solver != source] object type that gives mesh info and domain parameters",
  "geom": {
    "#mesh": "[required] object type that sets the mesh file type, name, and globally scaling factor",
    "mesh": {
      "type": "gmsh",
      "file": "sphere_1.msh",
      "scale": 2.5E-8
    },
		"#qd_tri": "[optional] string type that specifies the quadrature rule on triangular surface of all domains from: tri_gl[n] with n = 1|3|4|7|13, by default n=13",
		"qd_tri": "tri_gl13",
		"#qd_tetra": "[optional] string type that specifies the quadrature rule on tetra volume of all domains from: tetra_gl[n] with n = 1|4|11, by default n=4",
		"qd_tetra" "tetra_gl4",
    "#domain": "[required] array of object type that sets the parameters of the corresponding domain. Note that the 1st one is always used to defined the domain where the incident field is defined.",
    "domain": [
      {
        "#media": "[required] number (integer) type for the media id of which the physical properties are applied. Note that the media index starts at 0 in JSON array, i.e., .physics.media[0] is the first media.",
        "media": 0,
        "#surface": "[required] array of number type for the id[s] of the surface[s] enclose[s] the domain. The integers may be negative, which is interpreted so that the orientation of the corresponding surface is reversed. It is required that the boundary surface of a domain is oriented so that the normals point into the domain. Thus, if a surface id appears in one domain.surface expression, it must appear in another domain.surface with a minus sign.",
        "surface": [-11, -12],
        "#volume": "[required] array of number type for the volume id[s]",
        "volume": [],
        "qd_tri": "tri_gl13",
        "qd_tetra": "tetra_gl4"
      },
      {
        "media": 2,
        "surface": [11, 12],
        "volume": []
      }
    ]
  },
  "#physics": "[required] object type that specifies the physical parameters and further contains subnodes",
  "physics": {
    "#media": "[required] array with each element being an object that specifies the permittivity and permeability of the media",
    "media": [
      {
        "comment": "surrounding medium",
        "#name": "[required] string type (maximum length, 32 characters) that labels this media",
        "name": "air/vacuum",
        "#ri|epsilon": "[required] object type that reads permittivity by either refractive index or by epsilon. When both are present, only the ri node is actively used.",
        "ri": {
          "#method": "[required] string type that controls the way to get value: [value|table|model]",
          "method": "value",
          "#real": "[optional if method==value] number type that gives the real part, default value is 1.0",
          "real": 1,
          "#imag": "[optional if method==value] number type that gives the imaginary part, default value is 0.0",
          "imag": 0,
          "#ref_file": "[required if method==table] string type that gives the file name under [ref/] folder, i.e., [ref_file].ref",
          "ref_file": "agjc",
          "#model": "[required if method==model] string type that specifies the model type: [L4|D+2CP|Drude if name==gold|Au|silver|Ag]",
          "model": "gold"
        },
        "#mu": "[optional] object type that reads permeability by mu, default value is 1.0",
        "mu": {
          "method": "value",
          "real": 1
        }
      },
      {
        "comment": "scatter 1",
        "name": "gold",
        "epsilon": {
          "method": "model",
          "model":  "D+2CP"
        },
        "#nonlinear": "[optional] object type that specifies the nonlinearity",
        "nonlinear": {
          "#chi2": "[optional] object type that controls second harmonic generation, SHG",
          "chi2": {
            "#type": "[optional] string type that specifies [surface|volume] contribution",
            "type": "surface",
            "#value": "[optional] array type with each element being a number type for [\chi2_nnn, \chi2_ntt, \chi2_ttn]",
            "value": [1, 0, 0]
          }
        }
      }
    ],
    "#source": "[required] array with each element being an object type that specifies the input field",
    "source": [
      {
        "#type": "[required] string type that switches among source type: [focus|pw]",
        "type": "focus",
        "focus": {
          "paraxial": true,
          "input": "lg",
          "norm": true,
          "focal": 1e-3,
          "waist": 1e-3,
          "na": 0.8,
          "#E0":	10,
          "lg": {
            "pol": "r",
            "n": 0,
            "l": 3
          }
        }
      }
    ]
  },
  "#simulation": "[required] object type that controls the simulation method",
  "simulation": {
    "#solver": "[required] object type that sets solver's parameters",
    "solver": {
      "#name": "[required] string type that telss the solver's name: [source|scattered|mode]",
      "name": "source",
      "#parallel": "[optional] boolean type that enables parallel computation or not, default [false]",
      "parallel": true,
      "#memory": "[optional] boolean type that enables keeping some varaibles in memory for post-processing in solution part, default [false]",
      "memory": false
    },
    "#wavelength": "[optional] object type that takes wavelengths",
    "wavelength": {
      "#method": "[required] string type that sets the way to take wavelength values: [list|range|zrange]",
      "method": "list",
      "#value": "[required] array of number type that take parameters to generate wavelength values",
      "value": [1060e-9]
    },
    "#frequency": "[optional] object type that takes frequencies",
    "frequency": {
      "#method": "[required] string type that sets the way to take frequency values: [list|range|zrange]",
      "method": "list",
      "#value": "[required] array of number type that take parameters to generate frequency values",
      "value": [1060e-9]
    }
  },
  "solution": {
    "#base": "[optional] object type that controls the way that the basic variables are saved",
    "base": {
      "#save": "[optional] number (integer) type at which wavelength index the basic solution data are saved. 1: save at the first wavelength, [2, 3] save at the 2nd and 3rd wavelengths, [-2, -1] save at the 1st and 2nd wavelengths in the reversed order",
      "save": [1, -2]
    },
    "#focal": "[optional] object type that sets where to calculate the electric fields",
    "focal": [
      {
        "#label": "[required] string type that labels the name of the cut plane where electric fields are calculated",
        "label": "xy",
        "#isrc": "[required] number type that sets the source index used in simulation: [1|...|number of sources]",
        "isrc": 1,
        "#wl": "[required] number type that sets the simulation wavelength",
        "wl": 1060e-9,
        "#z": "[required] array of number type that sets z-coordinate: [npoints, z_start, z_end]",
        "z": [1,  0],
        "#x": "[required] array of number type that sets x-coordinate: [npoints, x_start, x_end]",
        "x": [151, -2.4e-6,  2.4e-6],
        "#y": "[required] array of number type that sets y-coordinate: [npoints, x_start, x_end]",
        "y": [153, -2.5e-6,  2.5e-6]
      }
    ]
  }
}
```

#### 2. Read parameters from commands

This section explains the commands implemented in the main interface of the program. Note that the usage of the program should not be limited to the use of this interface. The interface only provides actions for the most common use of the program as was required during the development phase. The user may add new functionality to the source code and if necessary, add a new interface command for this.

Wavelengths are used to fix the frequency of the time-harmonic problem. This is mostly because usually resonances are intuitively related to size of geometrical features. Thus the wavelengts referenced below always correspond to the wavelength in vacuum.

In the following, brackets/braces are used to denote formal arguments of commands. Do not include the brackets in the actual arguments. A bracketed argument is to be replaced by only single actual argument and a braced argument can be replaced by multiple actual arguments.

The implementation uses SI units, so all parameters should represent values in this unit system.

To write a comment, simply put the hash mark # before your desired comment
```
Example: # This is a comment
```

COMMAND: `name [str]`  
DESCRIPTION:  
Sets the name of the computation. The name is used mainly to create output filenames. Avoid spaces and special characters in the name as this might cause trouble depending on the filesystem in use.

Example: name 'myscatteringproblem'

COMMAND: `mesh [filename] [scale]`  
DESCRIPTION:  
Loads the surface mesh that is used to construct the basis functions. The `[filename]` should point to a valid mesh file of some of the following formats:

- Gmsh mesh (*.msh)
- Netgen neutral mesh format (*.nmf)

The `[scale]` is a real number that is used to scale the points in the meshfile into points in the computational domain, which has units of meters.

Example: `mesh 'sphere.msh' 1E-9`

COMMAND: `quad [tri] [tetra]`  
DESCRIPTION:  
Quadrature rule for triangle and tetrahedron elements. Possible choises for `[tri]` are

```
tri_gl1
tri_gl3
tri_gl4
tri_gl7
tri_gl13
```
and for `[tetra]`

```
tetra_gl1
tetra_gl4
tetra_gl11
```
The number refers to the number of points.

COMMAND: `wlrg [nwl] [first] [last]`  
DESCRIPTION:  
Allocates datastructures for amount of `[nwl]` computations, which will have wavelengths that are equally spaced from `[first]` to `[last]`, including the endpoints.

Example: `wlrg 10 400E-9 1200E-9`

COMMAND: `wlls [nwl] {wavelengths}`  
DESCRIPTION:  
Allocated datastrctures for amount of `[nwl]` computations, which will have wavelengths given by the space separated list `{wavelengths}`.

Example: `wlls 3 300E-9 1.7E-6 2.3E-6`

COMMAND: `nsrc [n]`  
DESCRIPTION:  
Allocates space for `[n]` sources. Must be called prior to ssrc.

COMMAND: `ssrc [index] [type] {parameters}`  
DESCRIPTION:  
Define the source of index `[index]` in the scattering problem. Currently following parameters are valid:

`[type] {parameters}` ; description

```
'pw' [theta] [phi] [psi] ; A time-harmonic linearly polarized plane-wave
'focus_rad' [focal] [waist] [na] [normalize] ; Focused radially polarised beam.
'focus_azimut' [focal] [waist] [na] [normalize] ; Focused azimutally polarized beam
'focus_x' [focal] [waist] [na] [normalize] ; x-polarised focused gaussian beam
'focus_y' [focal] [waist] [na] [normalize] ; y-polarised focused gaussian beam
'focus_hg01' [focal] [waist] [na] [normalize] ; x-polarised focused Hermite-gaussian beam, mode HG_01
'dipole' [x] [y] [z] [dx] [dy] [dz] ; An electric dipole source at point (x,y,z) with dipole moment (dx,dy,sz)
```

The angle theta denotes the elevation angle i.e. angle between the z-axis and direction of propagation k. Angle phi denotes the azimuthal angle measured with respect to x-axis towards the y-axis. Angle psi denotes the polarization angle, measured from `k x z`. Angles are in degrees. To be more precise:

```
kx = k0*sin(theta)*cos(phi)
ky = k0*sin(theta)*sin(phi)
kz = k0*cos(theta)
```
where `k0 = 2*pi/lambda`

```
Ex = E0*(cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi))
Ey = E0*(-cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi))
Ez = E0*sin(psi)*sin(theta)
```
Parameters `[focal]`, `[waist]` and `[na]` denote focal length, beam waist and numerical aperture, respectively. If parameter `[normalize]` is `'true'`, then the fields are scaled by the maximum electric field at the focus. If parameter `[normalize]` is `'false'`, then the fields are normalized to the incident plane-wave incident of the focusing lens.

COMMAND: `scan [sz] [d] [npt]`  
DESCRIPTION:  
Defines a special set of sources, where a source defined by a single `'ssrc'` command is scanned over a square area of side length `[d]` over plane `z = [sz]`. The number of scan points is `[npt]*[npt]`. Example usage:

```
nsrc 1
ssrc 1 'focus_rad' 1e-3 0.8e-3 0.8 false
scan 0 2e-6 20
```

COMMAND: `nmed [n]`  
DESCRIPTION:  
Allocates space for `[n]` media descriptors.

COMMAND: `smed [mindex] [type] [method] {parameters}`  
DESCRIPTION:  
Sets the medium of index `[mindex]`. The medium type `[type]` may be one of

```
linear
nlsurf
nlbulk_nonlocal
```

The paramter `[method]` describes how a property is set. Valid choices are

`value`
`file`

For type `'linear'`, the `{parameters}` descrive the complex index of refraction. If method `'value'` is used, then `{parameters}` is the numeric value of the complex index of refraction in Fortran syntax and the value is used for all wavelengths. If method is `'file'`, then the indices of refraction are read from a file, whose name is given as `{parameters}`. The file must contain three columns: wavelength in meters, real part of refractive index, imaginary part of refractive index (non-negative).

For type `'nlsurf'`, only method `'value'` is valid. Then `{parameters} = [nnn] [ntt] [ttn]` where the three numbers are the corresponding second-order surface susceptiblity components for isotropic surfaces.

For type `'nlbulk_nonlocal'`, only method `'value'` is valid. Then `{parameters} = [zzz]` which is the zzz-component of the dipolar bulk second-order susceptibility (other components are zero).

COMMAND: `ndom [n]`  
DESCRIPTION:  
Sets the number of domains to `[n]`.

COMMAND: `sdom [dindex] [nsurf] [surf_ids] [nvol] [vol_ids] [refind_file] [gf_index]`  
DESCRIPTION:  
Sets domain parameters for domain of index `[dindex]`. Command ndom must be issued before this one. `[nsurf]` is the number of surfaces whose union constitutes the boundary of the domain in the mesh file. `[surf_ids]` are integers corresponding to the physical id:s of the surfaces. The integers may be also negative, which is interpreted so that the orientation of the corresponding surface is reversed. `[nvol]` in the number of volumes whose union constitutes the domain. This may be zero if the volume information is not required. `[vol_ids]` are integers corresponding to the physical id:s of the volumes (orientation is arbitrary). Parameter `[refind_file]` is a file name of the refractive index data to be used for the domain. Parameter `[gf_index]` is an index to a Green's function that has been loaded by command ipgf. A value of `-1` corresponds to the nonperiodic Green's function.

It is required that the boundary surface of a domain is oriented so that the normals point into the domain. Thus, if a surface id appears in one sdom expression, it must appear in another one with a minus sign.

It is important to note that the domain, where the incident field is defined, always corresponds to `dindex = 1`.

COMMAND: `fdom`  
DESCRIPTION:  
Finalizes the domain setup by dividing the given mesh into proper submeshes while keeping book on the edge connectivity. It also orients the basis function on the submeshes properly. This should be called right before solv.

COMMAND: `sbnd [id] [bndname]`  
DESCRIPTION:  
Assigns boundary values for edges with physical group `[id]`. The last argument is a string from the following list

```
prdx1
prdx2
prdy1
prdy2
```

which denote the two ends of the unit cell in z = constant plane. E.g. prdx1 corresponds to edge x = -px/2, where px is the period in the direction of the x-axis. Basis coefficients for edges on prdx1 and prdy1 are removed from the system and added to coefficients related to edges prdx2 and prdy2 with proper Bloch phase shifts. These conditions must be set if the mesh touches the unit cell boundaries.
This routine must be called right after the 'mesh' command.

COMMAND: `symm [nsubgroups] {names}`  
DESCRIPTION:  
Defines the symmetry of the geometry. The symmetry is described by a group, which is generated from the given subgroups. The number of these subgroups is `[nsubgroups]`. A number of `[nsubgroups]` names must follow. The list of valid names is


id	identity
`mxp`	mirror symmetry with respect to the x-plane
`myp`	mirror symmetry with respect to the y-plane
`mzp`	mirror symmetry with respect to the z-plane
`r[n]`	[n]-fold rotation symmetry with respect to the z-axis

It is required that the symmetry group is commutative. Thus, for examples, command

`symm 2 mzp r3`

is valid but command

`symm 2 mxp r3`

is not valid. No error checking is done, so this is user's responsibility.

COMMAND: `solv`  
DESCRIPTION:  
Solves the problem.

COMMAND: `nfms [wlindex] [srcindex] [dindex]`  
DESCRIPTION:  
Computes the electric and magnetic fields corresponding to source `[srcindex]` on the surface of domain `[dindex]` of the particle at wavelength determined by `[wlindex]`. Produces a msh-file which can be inspected in gmsh.

COMMAND: `crst`  
DESCRIPTION:  
Computes scattering and absoprtion cross-sections at all available wavelengths and sources.

COMMAND: `rcst [wlindex] [srcindex] [ntheta] [nphi]`  
DESCRIPTION:  
Computes the bi-static radar cross-section (a.k.a. scattering power per unit solid angle) using the solution denoted by `[wlindex]` and source of index `[srcindex]`. Parameters `[ntheta]` and `[nphi]` denote how many angular points are evaluated for elevation and azimuthal angles respectively.

COMMAND: `rcs2 [wlindex]`  
DESCRIPTION:  
Computes the RCS integrated over the solid angle determined by a focused beam numerical aperture (in reflection). The solution at wavelength corresponding to `[wlindex]` is used.

COMMAND: `npgf [n]`  
DESCRIPTION:  
Allocate memory for `[n]` Green's functions. This must be called if periodic problems are to be solved.

COMMAND: `ipgf [index] [filename]`  
DESCRIPTION:  
Loads the precalculated coefficients for a periodic Green's function from file given by `[filename]` and assigns the data to positive integer `[index]`. Command npgf must be issued prior to calling ipgf.

The Green's function data is pre-computed by a MATLAB code made by the author. See the function file gpWlRange.m.

COMMAND: `diff [srcindex] [dindex] [orderx] [ordery] [pol] [polangle]`  
DESCRIPTION:  
Computes the diffracted power in the given order in a periodic problem. Integer `[srcindex]` refers to the excitation source. The integer `[dindex]` denotes the domain, where the calculation is done, i.e., which field expressions are used. Integers `[orderx]` and `[ordery]` denote the diffraction orders. `[pol]` may be '`true`' or '`false`' to specify whether a linear polarization filter is used when recording the diffracted wave. `[polangle]` is the polarizer angle in degrees (see `diffr.f90` for the actual definition of the angle).

It is assumed that the incident plane-wave propagates from half-space z>0 to half-space z<0.

##### Example input file

Say you want to compute the cross-sections of a sphere. Then pass in this input file:

```
name 'sphere'  
symm 1 id  
mesh 'unitsphere.msh' 200d-9  
wlrg 100 400d-9 1200d-9  
nmed 2  
smed 1 linear value (1.0, 0.0)  
smed 2 linear file aujc.ref  
ndom 2  
sdom 1 1 15 0 1 -1  
sdom 2 1 -15 0 2 -1  
fdom  
nsrc 1  
ssrc 1 'pw' 0 90 0  
solv  
crst  
exit  
```

Here the physical id of the surface is 15 and is originally oriented so that the normals point into the surrounding space.

A more involved case would be a periodic problem, where a particle stands on a substrate.
The input would look something like this:

```
name 'myproblem'  
symm 1 id  
wlrg 100 4e-7 10e-7  
mesh 'mymesh.msh' 1e-9  
nsrc 1  
ssrc 1 'pw' 180 0 0  
nmed 3  
smed 1 linear value (1.0, 0.0)  
smed 2 linear file aujc.ref  
smed 3 linear file silica.ref  
npgf 2  
ipgw 1 'vacuum.pgf'  
ipgw 2 'silica.pgf'  
ndom 3  
sdom 1 2 187 188 1 1  
sdom 2 2 -188 -189 2 -1  
sdom 3 2 -187 189 3 2  
fdom  
solv  
diff 1 3 0 0 'false' 0  
exit  
```

The glass substrate occupies the negative z half-space and the plane-wave is incident in the negative z direction.

##### File format specifications

Here brackets denote a row or a matrix of data and braces denote an arbitrary sequence of data.

EXTENSION: `crs`
DESCRIPTION:
Scattering and absorption cross-sections.
CONTENT:
`[wavelengths (float)] [scattering cross-sections (float)] [absorption cross-sections (float)]`

EXTENSION: `dif`
DESCRIPTION:
Diffracted intensity in zeroth order normalized to intensity of the incident plane-wave.
CONTENT:
`[wavelengths (float)] [intensity (float)]`
