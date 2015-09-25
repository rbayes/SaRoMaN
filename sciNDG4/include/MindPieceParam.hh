//MindPieceParam
//
//Define complete volume of mind according to parameterisation of Piece positions.

#ifndef MINDPIECEPARAM_HH
#define MINDPIECEPARAM_HH 1

#include <vector>

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;

class MindPieceParam : public G4VPVParameterisation
{
public:

  MindPieceParam(G4int npieces, G4double pieceLen);

  virtual ~MindPieceParam();

  virtual void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const;
  
private :
  
  std::vector<double> y;

};

#endif
