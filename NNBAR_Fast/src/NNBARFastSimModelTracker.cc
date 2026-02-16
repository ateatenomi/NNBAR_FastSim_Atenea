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
/// \file NNBARFastSimModelTracker.cc
/// \brief Implementation of the NNBARFastSimModelTracker class

#include "NNBARFastSimModelTracker.hh"
#include "NNBAREventInformation.hh"
#include "NNBARPrimaryParticleInformation.hh"
#include "NNBARSmearer.hh"
#include "NNBAROutput.hh"

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include "G4RegionStore.hh" //ate 3.2.26
#include "G4GeometryManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "Randomize.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4PathFinder.hh"
#include "G4FieldTrack.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NNBARFastSimModelTracker::NNBARFastSimModelTracker( G4String aModelName, 
  G4Region* aEnvelope, NNBARDetectorParametrisation::Parametrisation aType ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( aType ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NNBARFastSimModelTracker::NNBARFastSimModelTracker( G4String aModelName, 
                                                    G4Region* aEnvelope ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( NNBARDetectorParametrisation::eNNBAR ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NNBARFastSimModelTracker::NNBARFastSimModelTracker( G4String aModelName ) :
  G4VFastSimulationModel( aModelName ), fCalculateParametrisation(),
  fParametrisation( NNBARDetectorParametrisation::eNNBAR ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NNBARFastSimModelTracker::~NNBARFastSimModelTracker() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool NNBARFastSimModelTracker::IsApplicable( const G4ParticleDefinition& 
                                                                   aParticleType ) {
  return aParticleType.GetPDGCharge() != 0;  // Applicable for all charged particles
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool NNBARFastSimModelTracker::ModelTrigger( const G4FastTrack& /*aFastTrack*/ ) {
  return true;  // No kinematical restrictions to apply the parametrisation
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NNBARFastSimModelTracker::DoIt( const G4FastTrack& aFastTrack,
                                     G4FastStep& aFastStep ) {

  //G4cout << " ________Tracker model triggered _________" << G4endl;

  // Calculate the final position (at the outer boundary of the tracking detector)
  // of the particle with the momentum at the entrance of the tracking detector.

  //ate25 -  Mark acceptance
  auto track2 = aFastTrack.GetPrimaryTrack();

  //if (track2->GetParentID() == 0) {  // primary only
  if (true) {  // we also want muons (Ate) 29.1.26
    auto info = static_cast<NNBAREventInformation*>(
        G4EventManager::GetEventManager()->GetUserInformation());

    if (!info->GetPrimaryHitTracker()) {//Gets activated when the particle first enters the volume
        info->SetPrimaryHitTracker(); 
        info->SetEntryPos(track2->GetPosition());//Records entry position (Ate 3.2.26)
    }
  }



  G4int pdgID = 0;
  pdgID = aFastTrack.GetPrimaryTrack()-> GetDefinition()->GetPDGEncoding();

  if (( abs(pdgID) != 211 )and( abs(pdgID) != 2212 )){
  //if (( abs(pdgID) != 211 )and( abs(pdgID) != 2212 )and( abs(pdgID) != 13)) { //Ate 29.1 added muon
 // Kill the parameterised particle at the entrance of the electromagnetic calorimeter
  aFastStep.KillPrimaryTrack();
  aFastStep.ProposePrimaryTrackPathLength( 0.0 );
  }
  
  G4double KE = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();
  G4cout<<"Energy = "<<KE<<G4endl;


//Ate 29.1.26 - Duplicate
  G4Track track = * aFastTrack.GetPrimaryTrack();
  G4FieldTrack aFieldTrack( '0' );
  G4FieldTrack endTrack( 'a' );
  G4FieldTrackUpdator::Update( &aFieldTrack, &track );
  
  G4double Edep_pion = 0.0;

  //if ((abs(pdgID) == 13)){ //particle equal 13-muon ate 29.1.26
    
  //}

  if (( abs(pdgID) == 211)or( abs(pdgID) == 2212 )) {  // particle equal 211-pion+ (WVS) and 2212-proton (WVS) 2025-07-08
    ///*Ate 29.1.26
    G4Track track = * aFastTrack.GetPrimaryTrack();
    G4FieldTrack aFieldTrack( '0' );
    G4FieldTrackUpdator::Update( &aFieldTrack, &track );

    G4FieldTrack endTrack( 'a' );
    //*/ 
    G4double dEdx_rho = NNBARSmearer::Instance()->
         BetheBloch(aFastTrack.GetPrimaryTrack()-> GetDefinition(), 
                    KE);
    //std::cout << "dexrho: " << dEdx_rho << std::endl;


  G4double retSafety = -1.0;
  ELimited retStepLimited;
  G4double currentMinimumStep = 1.0*cm;  // Temporary: change that to sth connected
                                         // to particle momentum.
  G4PathFinder* fPathFinder = G4PathFinder::GetInstance();
  /*G4double lengthAlongCurve = */ 
  fPathFinder->ComputeStep( aFieldTrack,
                            currentMinimumStep,
                            0,
                            aFastTrack.GetPrimaryTrack()->GetCurrentStepNumber(),
                            retSafety,
                            retStepLimited,
                            endTrack,
                            aFastTrack.GetPrimaryTrack()->GetVolume() );

  // Place the particle at the tracking detector exit 
  // (at the place it would reach without the change of its momentum).
    //Calculate the energy depositied by a generic muon given the path length 
     G4ThreeVector posInitial = aFastTrack.GetPrimaryTrack()->GetPosition();
     G4ThreeVector posFinal = endTrack.GetPosition();

    //REGION CHECK (ate 3.2.26)
    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4VPhysicalVolume* finalVolume = navigator->LocateGlobalPointAndSetup(posFinal);
    G4Region* tpcRegion = G4RegionStore::GetInstance()->GetRegion("EM_tpc_region");
    if (finalVolume->GetLogicalVolume()->GetRegion() != tpcRegion) {
      // Particle left TPC
      auto info = static_cast<NNBAREventInformation*>(G4EventManager::GetEventManager()->GetUserInformation());
      info->SetExitPos(posFinal);
      G4double totPathLength = (posFinal-(info->GetEntryPos())).mag()/cm;
      info->SetTotPathLength(totPathLength);
      G4cout<<"DEBUG exit position: "<<totPathLength<<G4endl;
    }


     G4double pathLength = (posFinal - posInitial).mag()/cm;
     Edep_pion = ( dEdx_rho*pathLength )/MeV;

    /* G4cout << "DEBUG path=" << pathLength
       << " dEdx=" << dEdx_rho
       << " Edep=" << Edep_pion
       << G4endl; */
     /*if (Edep_pion <= 0.0 || !std::isfinite(Edep_pion)) { //ate jan25
  // No valid energy deposition → skip smearing & histogram
  aFastStep.ProposePrimaryTrackFinalPosition(
      aFastTrack.GetPrimaryTrack()->GetPosition()
  );
  aFastStep.ProposePrimaryTrackFinalKineticEnergy(KE);
  return;
}*/
    //DEBUG UNIT OF PATH LENGTH
    /*
    G4double myValue = pathLength;
  G4cout << "Raw value: " << myValue << G4endl;
  G4cout << "Value/mm: " << myValue/mm << G4endl;
  G4cout << "Value/cm: " << myValue/cm << G4endl;
  G4cout << "Value/MeV: " << myValue/MeV << G4endl;
  */


    NNBAROutput::Instance()->SetParticleTrace(pathLength); //WVS - 11/06/2025

    //Ate 28.1.25
    // Propagate kinetic energy, clamping to zero to avoid unphysical negative values  
    G4double KE_out = std::max(KE-Edep_pion,0.0 );

    // Define measured dE/dx from deposited energy and actual track length  
    double dEdx_meas= Edep_pion/pathLength;

    // Reject ultra-short tracks where dE/dx is unstable or meaningless  
    G4double L_min= 1*cm;
    //if (pathLength<L_min) return;

    aFastStep.ProposePrimaryTrackFinalPosition( endTrack.GetPosition() );
    aFastStep.ProposePrimaryTrackFinalKineticEnergy(KE_out);
    
    if ( !aFastTrack.GetPrimaryTrack()->GetParentID() ) {

          auto event =
        G4RunManager::GetRunManager()->GetCurrentEvent();

          auto info = static_cast<NNBAREventInformation*>(
        event->GetUserInformation());

          info->SetPrimaryHitTracker();
    }

    /*Ate 28.1.25
    //if (KE < Edep_pion) {
    if (false) {
        std::cout<<"Killed"<<std::endl;
        aFastStep.KillPrimaryTrack();
        //aFastStep.ProposePrimaryTrackPathLength( 0.0 );
    }
    
    else {
        G4double DeltaKE = KE - Edep_pion;
        //KE = Edep_pion; //nonsense
        //std::cout << ">>>>>>>> Tracker DeltaKE (MeV) " << ": " <<  DeltaKE << std::endl;
        //Place the particle at the detector exit 
        //(at the place it would reach without the change of its momentum).
        aFastStep.ProposePrimaryTrackFinalPosition( endTrack.GetPosition() );
        aFastStep.ProposePrimaryTrackFinalKineticEnergy(DeltaKE);

        if ( !aFastTrack.GetPrimaryTrack()->GetParentID() ) {

          auto event =
        G4RunManager::GetRunManager()->GetCurrentEvent();

          auto info = static_cast<NNBAREventInformation*>(
        event->GetUserInformation());

          info->SetPrimaryHitTracker();
        }
    }*/







    }
  // Consider only primary tracks (do nothing else for secondary charged particles)





  G4ThreeVector Porg = aFastTrack.GetPrimaryTrack()->GetMomentum();
  if ( ! aFastTrack.GetPrimaryTrack()->GetParentID() ) {
    

    G4ThreeVector Pos = aFastTrack.GetPrimaryTrack()->GetPosition();
    //std::cout<<"Pos x"<< Pos.x()<<std::endl;
    G4double time = aFastTrack.GetPrimaryTrack()->GetGlobalTime();

    NNBAREventInformation* info = (NNBAREventInformation*) 
                            G4EventManager::GetEventManager()->GetUserInformation();


    
    if ( info->GetDoSmearing() ) {
      //std::cout<<"I am in GetDoSmearing"<<std::endl;
      // Smearing according to the tracking detector resolution
      G4double res = fCalculateParametrisation->
        GetResolution( NNBARDetectorParametrisation::eTRACKER, 
                       fParametrisation, KE, pdgID );
      G4double med = fCalculateParametrisation->
        GetMedian( NNBARDetectorParametrisation::eTRACKER, 
                   fParametrisation, KE, pdgID );
      G4double eff = fCalculateParametrisation->
        GetEfficiency( NNBARDetectorParametrisation::eTRACKER,
                       fParametrisation, KE );
      //G4ThreeVector Psm;
      //Psm = NNBARSmearer::Instance()->
      //                           SmearMomentum( aFastTrack.GetPrimaryTrack(), res );
      //NNBAROutput::Instance()->FillHistogram( 0, ((Psm.mag()/MeV) / (Porg.mag()/MeV)) );
      //std::cout << "prueba " << ": " <<  Psm.mag()/MeV  << std::endl;

      //G4cout << "DEBUG: Edep_pion = " << Edep_pion/MeV << G4endl;
      G4double Psm = NNBARSmearer::Instance()->
                                 SmearEnergy( aFastTrack.GetPrimaryTrack(), res, med, Edep_pion );
    NNBAROutput::Instance()->FillHistogram( 0, ((Psm/MeV) / (Edep_pion/MeV))); // ((Psm.mag()/MeV) / (Porg.mag()/MeV)) );

    NNBAROutput::Instance()->SaveTrack( NNBAROutput::eSaveTracker,
                                         0,
                                         pdgID,
                                         Edep_pion/MeV,
                                         Pos/mm,
                                         res,
                                         eff,
                                         Psm/MeV,                                         
                                         time/ns /*,
                                         Volume*/);  //Volume - WVS - 14/05/2025

      //to record only one value for energy smearing on total trace - WVS- 24/06/2025
      if (!NNBAROutput::Instance()->IsSmearingSelect()) {
        G4double SmearingGaussValue = NNBARSmearer::Instance()->SmearValue( res );
        NNBAROutput::Instance()->SetSmearingValueSelect(SmearingGaussValue);
      }

      // Setting the values of Psm, res and eff
      /*( (NNBARPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerMomentum( Psm );
      ( (NNBARPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerResolution( res );
      ( (NNBARPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerEfficiency( eff );*/
    G4ThreeVector posInitial = aFastTrack.GetPrimaryTrack()->GetPosition();
    G4ThreeVector posFinal = endTrack.GetPosition();
    G4double pathLength = (posFinal - posInitial).mag()/cm;
      
    } else {
      // No smearing: simply setting the value of Porg
      ( (NNBARPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerMomentum( Porg );
        aFastStep.ProposeTotalEnergyDeposited( KE );
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

