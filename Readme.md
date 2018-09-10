# Sound propagation in sonic crystals
The following routines solve different Helmholtz problems for the acoustic pressure. It is meant to analyze acoustic waves in periodic media.

All the routines use python3.
To run the 2D codes, you need: 
- The FEniCS library ([https://fenicsproject.org/](https://fenicsproject.org/))
- gmsh ([http://gmsh.info/](http://gmsh.info/))
- usual python packages. 


To run the 3D codes you do not need FEniCS but you need the bempp library([https://bempp.com/](https://bempp.com/)).

My set up is Ubuntu 16.04 with:
- fenics (2017.2.0)
- matplotlib (2.2.2)
- numpy (1.14.5)
- scipy (1.1.0)
- setuptools (39.2.0)
- bempp (3.3.3)

Here is a quick overview of the routines with the corresponding [folder/file]:

**1D [1D]:**
- Illustration of the dispersion relation in free space [/DispRel.py]
- Analytic dispersion relation [/Multilayer/multilayer.py]
- Transfer matrix method [/Multilayer/multilayer_TMM.py]
- Discrete periodic mass-spring system and Helmholtz resonator [/Helmholtz_res/1DPeriodic_HelmholtzResModel.py]
- If you have Jupyter, you can run the interactive notebook [/Helmholtz_res/HelmholtzResModel.ipynb]


**2D [2D-fenics]: FEM set up for eigenvalues and transmission problems**
- Brillouin zones in 2D for square and hexagonal periodicity [/Brillouin_zones.py]
- Infinite periodic medium analysis [/BG/FEM_2D.py]
- 1-direction transmission [/Transmission/SC_1Dir/FEM_2D_Trans.py]
- 2-directions transmission [/Transmission/SC_PML/FEM_PML.py]

**3D [3D-bempp]: BEM set up for sound hard scattering**
- Unit sphere scattering over several frequencies [/Sphere/sphere_loop.py]
- Sonic crystal scattering [/SC_Jasa/SC_bempp.py], see reference for the geometry: *Karimi, M., Croaker, P., & Kessissoglou, N. (2017). Acoustic scattering for 3D multi-directional periodic structures using the boundary element method. _The Journal of the Acoustical Society of America_, _141_(1), 313-323.*

Default cases are used if you run the codes as they are.
