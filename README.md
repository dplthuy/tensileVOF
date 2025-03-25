# tensileVOF
Code repository of the Tensile Force method as described in :

* D.P.L. Thuy, N.G. Deen, J.J.C. Remmers, G. Finotello, "Volume-of-fluid simulations of multiphase flows with high surface tension and curvature: a tensile force approach with pressure jump correction" (Under review)*

## Installation instructions
- Source code and test cases are intended for use with a compiled version of [OpenFOAM v2406](https://www.openfoam.com/news/main-news/openfoam-v2406). Compatibility with any other version of OpenFOAM is not guaranteed.
- Clone the repository to the `$WM_PROJECT_USER_DIR` of your OpenFOAM installation.
  ```
  cd $WM_PROJECT_USER_DIR
  git clone https://github.com/dplthuy/tensileVOF.git
  ```
- Compile the source code by executing `src/Allwmake` or by using `wmake` in the respective directories.
