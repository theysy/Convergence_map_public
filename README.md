# Information of developer
    1. Name: Seong-Yong Yoon
    2. Affiliation: Pohang university of science and technology (POSTECH)
    3. Advisor: Frederic Barlat
    4. E-mail: theysy@postech.ac.kr

# Brief description
    - The convergence map (CMAP) is a numerical tool to approximatley estimate the numerical stability and computational efficiency.
    - CMAP records the number of iteration necessary for convergence of an arbitrary stress update algorithm.
    - The darkness of color indicates the number of iteration.
    - Individual iteration data are recored on deviatroic plane.
    - CMAP_U2 is calculated based on MML_U2 meanwhile CMAP_U3 is based on MML_U3.

# How to use CMAP code
    1. Fill out the user material property file in UMAT_PROPS folder. Ex) PROPS_AA6022_YLD2K_HAH20.CSV
    2. Run 'CMAPS.for' using intel fortran or gfortran
    3. 'CAMP.csv' contrains the data for convergence map.
    4. Run 'Convg_Mapping.py' using python 3.6.

# Prerequisites
  1. Intel or GNU fortran
  2. Python 3.6
  3. Python packages: numpy, matplotlib

# Screeshots
1. Anisotropic yield function: **Yld2000-2d**
<p align="center"><img src="https://github.com/theysy/Convergence_MAP/blob/main/Screenshot/AA6022_YLD2K_CMAP.png"></p>

2. Anisotropic hardening: **HAH<sub>20</sub>**
<p align="center"><img src="https://github.com/theysy/Convergence_MAP/blob/main/Screenshot/AA6022_YLD2K_HAH20_CMAP.png"></p>
