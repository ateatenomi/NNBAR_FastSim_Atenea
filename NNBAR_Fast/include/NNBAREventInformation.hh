//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file NNBAREventInformation.hh
/// \brief Definition of the NNBAREventInformation class

#ifndef NNBAR_EVENT_INFORMATION_H
#define NNBAR_EVENT_INFORMATION_H

#include "G4VUserEventInformation.hh"
#include "globals.hh"
#include "G4ThreeVector.hh" 

/// Event information.
///
/// Describes the information that can be associated with a G4Event class object.
/// @author Anna Zaborowska

class NNBAREventInformation : public G4VUserEventInformation {
  public:
    
    /// A default constructor. Sets flag fDoSmearing to true.
    NNBAREventInformation();
    
    /// A constructor.
    /// @param aSmear The flag indicating if smearing should be done.
    NNBAREventInformation( G4bool aSmear );

    virtual ~NNBAREventInformation();
    
    /// Prints event information.
    virtual void Print() const;
    
    /// Sets the flag indicating if smearing should be done.
    /// @param aSmear A boolean flag.
    void SetDoSmearing( G4bool aSmear );
    
    /// Gets the flag indicating if smearing should be done.
    G4bool GetDoSmearing();

    //ate dec25: Add an acceptance flag
    void Reset() {
        fPrimaryHitTracker = false;
    }
    void SetPrimaryHitTracker() { fPrimaryHitTracker = true; }
    G4bool GetPrimaryHitTracker() const { return fPrimaryHitTracker; }

    void SetEntryPos(G4ThreeVector pos) {fEntryPos = pos;}
    G4ThreeVector GetEntryPos() const {return fEntryPos;}

    void SetExitPos(G4ThreeVector pos2) {fExitPos = pos2;}
    G4ThreeVector GetExitPos() const {return fExitPos;}

    // Generated kinematics
    void SetGenMomentum(G4double p) { fGenMomentum = p; }
    void SetGenTheta(G4double t) { fGenTheta = t; }
    void SetGenKE(G4double KE) {fGenKE = KE;}
    void SetGenPrimaryInitialX(G4double posx) {fPrimaryInitialX =posx;}
    void SetTotPathLength(G4double path) {fTotPathLength=path;}

    G4double GetGenMomentum() const { return fGenMomentum; }
    G4double GetGenTheta() const { return fGenTheta; }
    G4double GetGenKE() const {return fGenKE;}
    G4double GetGenPrimaryInitialX() const { return fPrimaryInitialX; }
    G4double GetTotPathLength() const { return fTotPathLength;}

  private:
    
    /// A flag indicating if smearing should be performed. 
    /// It is read by implementations of G4VFastSimulationModel.
    G4bool fDoSmearing;

    //ate dec25:
    G4bool fPrimaryHitTracker = false;
    G4double fGenMomentum = 0.;
    G4double fGenTheta    = 0.;
    G4double fGenKE = 0.;
    G4double fPrimaryInitialX=0.;
    G4ThreeVector fEntryPos;
    G4ThreeVector fExitPos;
    G4double fTotPathLength = 0.; 
};

#endif

