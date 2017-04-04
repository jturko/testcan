
#ifndef Hit_h
#define Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class Hit : public G4VHit
{
    
public:
    Hit();
    ~Hit();
    
    inline G4double GetEnergy() {return energy;}
    inline void SetEnergy(G4double enr) { energy = enr;}
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
private:
    G4double energy;
    
};

typedef G4THitsCollection<Hit> MyHitsCollection;

extern G4ThreadLocal G4Allocator<Hit>* MyHitAllocator;

inline void* B5EmCalorimeterHit::operator new(size_t)
{
    if (!MyHitAllocator)
        MyHitAllocator = new G4Allocator<Hit>;
    return (void*)MyHitAllocator->MallocSingle();
}

inline void B5EmCalorimeterHit::operator delete(void* aHit)
{
    MyHitAllocator->FreeSingle((Hit*) aHit);
}

#endif
