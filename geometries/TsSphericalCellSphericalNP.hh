//
// ********************************************************************
// *                                                                  *
// * This file is based on the TsSphericalCell example                *
// * from the TOPAS-nBio extensions to the TOPAS Simulation Toolkit.  *
// * The TOPAS-nBio extensions are freely available under the license *
// * agreement set forth at: https://topas-nbio.readthedocs.io/       *
// *                                                                  *              
// *  Modifications by Marc B. Hahn (2020)                            *
// *  Please report bugs to hahn@physik.fu-berlin.de                 *
// *  or on https://github.com/BAMresearch                            *    
// ********************************************************************
//
// A simple spherical cell with nanoparticles can be generated in a fast manner.
// For elliptical subcomponents use TsSphericalCellNP which is more flexible but slower.
// User has the option of including organelles: nucleus, mitochondria, cell membrane and/or nanoparticles.


#ifndef TsSphericalCellSphericalNP_hh
#define TsSphericalCellSphericalNP_hh

#include "TsVGeometryComponent.hh"
#include <vector>


class TsSphericalCellSphericalNP : public TsVGeometryComponent
{    
public:
	TsSphericalCellSphericalNP(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsSphericalCellSphericalNP();
	
	G4VPhysicalVolume* Construct(); 
    void ResolveParameters();
    
private:
    
    G4double CellRadius;
    G4double NucleusRadius;
    G4double MitoRadius;
    G4double MembraneThickness;
    
    G4RotationMatrix* rotationMatrix;
    G4VPhysicalVolume* pNucleus;
    
    std::vector<std::vector<G4double> >  CellCoordinates;
    std::vector<std::vector<G4double> >  TargetSphereCoordinates;

    std::vector<G4double>  tmpCoordinates;
    
    G4bool CheckOverlapOfSphereWithGeometryComponents(std::vector<std::vector<G4double> >& Coordinates, G4double r, G4double x, G4double y, G4double z);
    
    G4ThreeVector* AddSphereToCell(G4double radius);
    void AddCoordinates(std::vector<std::vector<G4double> >& Coordinates, G4double r, G4double x, G4double y, G4double z);   
    
    
};

#endif
