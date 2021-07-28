#ifndef FluxHit_h
#define FluxHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"

//---------------------------------------------------------------------------

class FluxHit : public G4VHit
{
public:
  
  FluxHit();
  ~FluxHit();
  FluxHit(const FluxHit&);
  const FluxHit& operator=(const
						FluxHit&);
  int operator==(const FluxHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  virtual void Draw();
  void Print();

protected:

  G4int                  fID; 
  G4ParticleDefinition*  fPDef; 
  G4double               fTime; 
  G4double               fEdep;
  G4ThreeVector          fMom;   
  G4ThreeVector          fPosPre; 
  G4ThreeVector          fPosPost; 
  G4int                  fAtomicMass;
  G4int                  fAtomicNumber;
  G4int                  fTrackID;
  G4int                  fParentID;

public:
  
  inline void SetID       (G4int i)                  { fID  = i;   };
  inline void SetPDef     (G4ParticleDefinition* pd) { fPDef = pd;  };
  inline void SetTime     (G4double t)               { fTime = t;   };
  inline void SetMomentum (G4ThreeVector m)          { fMom = m;    };
  inline void SetPrePosition (G4ThreeVector pos)     { fPosPre=pos;    };
  inline void SetPostPosition (G4ThreeVector pos)    { fPosPost=pos;    };
  inline void SetEnergy   (G4double de)              { fEdep = de; };
  inline void SetAtomicMass (G4int mm)               { fAtomicMass = mm; };//SetAtomicMass (for ions)
  inline void SetAtomicNumber (G4int nn)             { fAtomicNumber = nn; };//SetAtomicNumber 
  inline void SetTrackID    (G4int tt)               { fTrackID = tt; };//SetTrackID (number associated with each track
  inline void SetParentID   (G4int pp)               { fParentID = pp; };//SetParentID(numberassociated with eachprimaryparticle

  inline G4int                 GetID()               { return fID;  };
  inline G4int                 GetTrackID()          { return fTrackID;  };  
  inline G4int                 GetParentID()         { return fParentID;  };
  inline G4ParticleDefinition* GetPDef()             { return fPDef; };
  inline G4double              GetTime()             { return fTime; };
  inline G4ThreeVector         GetMom()              { return fMom; };
  inline G4ThreeVector         GetPosPre()           { return fPosPre; };
  inline G4ThreeVector         GetPosPost()          { return fPosPost; };
  inline G4double              GetEdep()             { return fEdep; };  
  inline G4int                 GetAtomicNumber()     { return fAtomicNumber; };  
  inline G4int                 GetAtomicMass()       { return fAtomicMass; };  
};

//---------------------------------------------------------------------------

typedef G4THitsCollection<FluxHit> FluxHitsCollection;

extern G4Allocator<FluxHit> FluxHitAllocator;


inline void* FluxHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) FluxHitAllocator.MallocSingle();
  return aHit;
}


inline void FluxHit::operator delete(void* aHit)
{
  FluxHitAllocator.FreeSingle((FluxHit*) aHit);
}

#endif

//---------------------------------------------------------------------------










