# m<img src="./img/magpie-dark.svg#gh-dark-mode-only" style="height:1ch;"/><img src="./img/magpie-light.svg#gh-light-mode-only" style="height:1ch;"/>gpie  C++

This is an experimental version of the magpie Python and MATLAB libraries translated for C++

The library uses the Spectra and Eigen C++ libraries which are included as submodules.

Unlike MATLAB [`eigs`](https://www.mathworks.com/help/matlab/ref/eigs.html) or the [SciPy equivalent](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html) the Spectra library _does not_ use the ARPACK fortran library. The ability to calculate a large number of eigenvalues and veecors is limited but the is not necessary at ths point for the purposes of the magpie library.


TODO: [Web Assembly build with emscripten](https://observablehq.com/@rreusser/eigen) 