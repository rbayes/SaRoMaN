//Parameterisation of mind piece positions.
#include "MindPieceParam.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

MindPieceParam::MindPieceParam(G4int npieces, G4double pieceLen)
{
  y.resize( npieces );
  
  for ( unsigned int i = 0; i < npieces; ++i )
    y[i] = ( (double)i - (double)npieces / 2.0 ) * pieceLen + pieceLen / 2.0;
  
}

MindPieceParam::~MindPieceParam()
{
}

void	MindPieceParam::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  physVol->SetTranslation( G4ThreeVector( 0.0, 0.0, y[ copyNo ] ) );
}
