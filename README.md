# JCARFAC
Julia wrappers for the [ccarfac](https://github.com/analogouscircuit/ccarfac) library
(an implementation of Richard Lyon's CARFAC model of the auditory periphery).

This module requires libcarfac.so, which can be compiled from the [ccarfac](https://github.com/analogouscircuit/ccarfac) repo. The ccarfac repository is included as a submodule of this project, so, after cloning JCARFAC repo, it is only necessary to execute
`git submodule init` and then `git submodule update`.  There is a makefile in the ccarfac project for compiling
libcarfac.  Note that `ccall`, the function for executing C code within Julia, requires a literal expression
for the location of the library being called.  You may have to modify the library location at each instance of `ccall`
in the Julia module. (If the file structure used in this repo is preserved and you are running a *nix-like
OS, it is likely that modifications will not be necessary.)

Once setup is complete, usage is very simple.  For a signal, `sig` (an array of 64-bit floats), and a sampling rate, `fs`, 
just call `JCARFAC.calcnap(sig, fs)`.  This will return a tuple `(nap, cfs)`, where `nap` will be a 2D array containing
the neural activity pattern, and `cfs` will contain the center frequencies for each channel.  Plotting via `heatmap(nap)`
will give you quick visual feedback. It is possible to modify the parameters of the model by declaring parameter structures
for the basilar membrane (`BMParams`), inner hair cells (`IHCParams`), and outer hair cells (`OHCParams`).  With these in
hand, simply call `calcnap(bmp, ihcp, ohcp, signal, fs)`, where `bmp`, `ihcp`, and `ohcp` are your parameter structures.
