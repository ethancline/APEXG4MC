#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "EnergyDepositHit.hh"

//---------------------------------------------------------------------------

G4Allocator<EnergyDepositHit> EnergyDepositHitAllocator;

EnergyDepositHit::EnergyDepositHit()
{
  fID   = 0;
  fEdep = 0;
  fTime = 0;
}

//---------------------------------------------------------------------------

EnergyDepositHit::~EnergyDepositHit()
{
}

//---------------------------------------------------------------------------

EnergyDepositHit::EnergyDepositHit(const EnergyDepositHit& right)
  :G4VHit()
{
  fID   = right.fID;
  fEdep = right.fEdep;
  fTime = right.fTime;

  fPosPre  + right.fPosPre;
  fPosPost  + right.fPosPost;
}

//---------------------------------------------------------------------------

const EnergyDepositHit& EnergyDepositHit::operator=(const EnergyDepositHit& right)
{

  fID   =right.fID;
  fEdep =right.fEdep;
  fTime = right.fTime;

  fPosPre  +right.fPosPre;
  fPosPost  +right.fPosPost;
  return *this;
}

//---------------------------------------------------------------------------

int EnergyDepositHit::operator==(const EnergyDepositHit&) const
{return 0;}

//---------------------------------------------------------------------------

void EnergyDepositHit::Draw()
{;}

//---------------------------------------------------------------------------

void EnergyDepositHit::Print()
{;}

//---------------------------------------------------------------------------











