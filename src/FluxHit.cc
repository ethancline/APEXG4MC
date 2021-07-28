#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "FluxHit.hh"

//---------------------------------------------------------------------------

G4Allocator<FluxHit> FluxHitAllocator;

FluxHit::FluxHit()
{
  fID   = 0;
  fPDef = 0;
  fTime = 0;
  fEdep = 0;
  fAtomicMass = 0;
  fAtomicNumber = 0;
  fTrackID = 0;   
  fParentID = 0; 
}

//---------------------------------------------------------------------------

FluxHit::~FluxHit()
{
}

//---------------------------------------------------------------------------

FluxHit::FluxHit(const FluxHit& right)
  :G4VHit()
{
  fID   = right.fID;
  fEdep = right.fEdep;
  fPDef = right.fPDef;
  fTime = right.fTime;
  fAtomicMass   = right.fAtomicMass;
  fAtomicNumber = right.fAtomicNumber;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;

  fMom  + right.fMom;
  fPosPre  + right.fPosPre;
  fPosPost  + right.fPosPost;
}

//---------------------------------------------------------------------------

const FluxHit& FluxHit::operator=(const FluxHit& right)
{

  fID   =right.fID;
  fEdep =right.fEdep;
  fPDef =right.fPDef;
  fTime =right.fTime;
  fAtomicMass   = right.fAtomicMass;
  fAtomicNumber = right.fAtomicNumber;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;

  fMom  +right.fMom;
  fPosPre  +right.fPosPre;
  fPosPost  +right.fPosPost;
  return *this;
}

//---------------------------------------------------------------------------

int FluxHit::operator==(const FluxHit&) const
{return 0;}

//---------------------------------------------------------------------------

void FluxHit::Draw()
{;}

//---------------------------------------------------------------------------

void FluxHit::Print()
{;}

//---------------------------------------------------------------------------











