## Notes

A work in progress to create a modern Fortran version of the filterSD optimization method.

## Original Fortran 77 Documentation

**FilterSD** is a package of Fortran 77 subroutines for solving nonlinear programming problems and linearly constrained problems in continuous optimization.

Please see the `README.pdf` file for instructions on compiling the source code.

This package does not require any Third Party software.

### README.pdf

The user needs to choose whether to use a sparse matrix or dense matrix data structure. To solve an NLP using a sparse matrix data structure, the subroutines in the following files are required
```
filterSD.f, checkd.f, glcpd.f, l1sold.f, shared.f, schurQR.f, sparseA.f, util.f
```
together with a user supplied driver program.

The file `schurQR.f` implements a Schur complement scheme for sparse matrix updates. This replaces an a previous file `sparseL.f` implementing Fletcher-Matthews updates, which is also included in the distribution. These files are interchangeable.

To solve an NLP using a dense matrix data structure, the subroutines in the following files are required
```
filterSD.f, checkd.f, glcpd.f, l1sold.f, shared.f, denseL.f, denseA.f, util.f
```
together with a user supplied driver program.

To solve an LCP using a sparse matrix data structure, the subroutines in the following files are required
```
glcpd.f, checkg.f, shared.f, schurQR.f, sparseA.f, util.f
```
together with a user supplied driver program.

To solve an LCP using a dense matrix data structure, the subroutines in the following files are required
```
glcpd.f, checkg.f, shared.f, denseL.f, denseA.f, util.f
```
together with a user supplied driver program.

Information on how to set up the driver program is contained in the files `filterSD.pdf` and `glcpd.pdf`. Examples of driver programs are provided in the files `hs106.f`, `hs106d.f`, `hs72.f` and `hs72d.f`. To solve a QP or LP, replace `glcpd.f` by `qlcpd.f` in the above. Usage of `qlcpd.f` is described at the head of the file and is similar to that for `glcpd.f`. To facilitate access to CUTEr NLP test problems, a driver program `driver.f` and associated subroutines in the file `user.f` is provided.

## Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/filterSD/). This was generated from the source code using [FORD](https://github.com/cmacmackin/ford).



