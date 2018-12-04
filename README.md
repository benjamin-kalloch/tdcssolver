# TDCSSolver

A solver application implementing the OpenFOAM API to compute a) the electrical field strength (E), b) the electrical current density (J) and c) electrical potential (ElPot) of a test case simulating transcranial direct current stimulation (tDCS). The application supports an arbitrary number of electrodes and conductivity values defined either as scalar values or tensor values. The application is described in the paper **A flexible open-source pipeline for simulating transcranial electric stimulation** by *Benjamin Kalloch, Pierre-Louis Bazin, Arno Villringer, Bernhard Sehm, and Mario Hlawitschka*.

## How to compile
The application compiles using the OpenFOAM v6 API:
1) You will need the source-pack release of OpenFOAM to be able to compile this solver application. You can download this release from [2].
2) Follow the instructions in [2] to compile and set up OpenFOAM.
3) If OpenFOAM is correctly set up, the application may be compiled simply by issuing the "wmake" command in the base directory of this repository. 

##### Useful links
[1] The OpenFOAM project: https://openfoam.org/

[2] Source pack release of OF: https://openfoam.org/download/6-source/
