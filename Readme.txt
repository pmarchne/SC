The following routines solve different Helmholtz problems for the acoustic pressure. It is meant to analyze acoustic waves in periodic media.

All the routines use python3.
To run the 2D codes, you need: 
- The FEniCS library (https://fenicsproject.org/)
- gmsh (http://gmsh.info/)
- usual python packages. 
To run the 3D codes you do not need FEniCS but you need the bempp library(https://bempp.com/).

My set up is Ubuntu 16.04 with:
- fenics (2017.2.0)
- matplotlib (2.2.2)
- numpy (1.14.5)
- scipy (1.1.0)
- setuptools (39.2.0)
- bempp (3.3.3)

Here is a quick overview of the routines [with the corresponding folder/file]:
1D [1D/]:
- Illustration of the dispersion relation in free space [DispRel.py]
- Analytic dispersion relation [Multilayer/multilayer.py]
- Transfer matrix method [Multilayer/multilayer_TMM.py]
- discrete periodic mass-spring system and Helmholtz resonator [Helmholtz_res/1DPeriodic_HelmholtzResModel.py]
- If you have Jupyter, you can run the interactive notebook [Helmholtz_res/HelmholtzResModel.ipynb]
###############################################################################################
2D [2D-fenics/]: FEM set up for different cases
- Brillouin zones in 2D for square and hexagonal periodicity [Brillouin_zones.py]
- Infinite periodic medium analysis [BG/FEM_2D.py]
- 1 direction transmission problem [Transmission/SC_1Dir/FEM_2D_Trans.py]
- 2 directions transmission problem [Transmission/SC_PML/FEM_PML.py]
###############################################################################################
3D [3D-bempp/]: BEM set up for differents cases
- sound-hard scattered sphere over several frequencies [Sphere/sphere_loop.py]
- sonic crystal scattering (from JASA 2017) [SC_Jasa/SC_bempp.py]
###############################################################################################

Default cases are used if you run the codes as they are.