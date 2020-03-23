# TOPAS-CellModels
1) Description:

Cell Models for TOPAS/Geant4 and the inclusion of nano particles in particle scattering simulations.

The C++ classes in this repository extend the functionality of the TOPAS (http://www.topasmc.org/) MC program, which is itself a wrapper of the Geant4 MCS Toolkit (http://geant4.org).


2) Installation:

Installation:

Navigate to the TOPAS extension directory:

  cd ~/topas_extensions/

Clone or download the sourcecode into your TOPAS extension directory:
 
  git clone https://github.com/MarcBHahn/TOPAS-CellModels.git
 
Change to your topas directory:

  cd ~/Topas/

Install it:
  
  cmake ./ -DTOPAS_EXTENSIONS_DIR=~/topas_extensions/TOPAS-CellModels &&  make -j4  


3) Usage:
