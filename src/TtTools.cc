// -*- C++ -*-
//
// Package:    TtTools
// Class:      TtTools
// 
/**\class TtTools TtTools.cc PhysicsTools/TtAnalysis/src/TtTools.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
//


// system include files
#include <memory>

// user include files
#include "TtTools.h"

#include "FWCore/Framework/interface/MakerMacros.h"

// constants, enums and typedefs

// static data member definitions


// constructors and destructor
using namespace edm;
using namespace std;
TtTools::TtTools()
{
   //now do what ever initialization is needed


  // get the associator parameters

}


TtTools::~TtTools()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtTools Analysis] Destructor called" << endl;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------


FreeTrajectoryState TtTools::getFTS(GlobalPoint GP, GlobalVector GV, int charge, 
                                  const AlgebraicSymMatrix66& cov, const MagneticField* field){


  GlobalTrajectoryParameters tPars(GP, GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

double TtTools::getEta(double vx, double vy, double vz ) {

      double va = sqrt( vx*vx + vy*vy + vz*vz );

      double theta = acos( vz/va );
      double eta = (-1.0)*log( tan(theta/2.0) )  ;
      return eta;
}

double TtTools::getdPhi( LorentzVector v1, LorentzVector v2 ) {

  double ab = (v1.Px()*v2.Px()) + (v1.Py()*v2.Py()) ;
  double al = sqrt( v1.Px()*v1.Px() +  v1.Py()*v1.Py() );
  double bl = sqrt( v2.Px()*v2.Px() +  v2.Py()*v2.Py() );
  double cosA = ab/(al*bl) ;
  double df = acos(cosA) ;

  return df;
}

double TtTools::getdR( LorentzVector v1, LorentzVector v2 ) {

  double df = getdPhi(v1, v2);

  double ah = getEta( v1.Px(), v1.Py(), v1.Pz() );
  double bh = getEta( v2.Px(), v2.Py(), v2.Pz() );
  double dh = ah - bh ;

  double dR = sqrt( (dh*dh) + (df*df) );

  return dR;
}

double TtTools::getdRy( LorentzVector v1, LorentzVector v2 ) {

  double df = getdPhi(v1, v2);

  double aY = getY( v1 );
  double bY = getY( v2 );

  double dY = aY - bY ;

  double dR = sqrt( (dY*dY) + (df*df) );

  return dR;
}

double TtTools::getY( LorentzVector v1 ){

    double ep = (v1.E() + v1.Pz()) / (v1.E() - v1.Pz()) ;
    double Y = -0.5*log( ep  ) ;
    return Y;
}

double TtTools::getInvMass( std::vector<LorentzVector> vlist ) {

    double cb[4] = {0.0};
    for ( std::vector<LorentzVector>::const_iterator it = vlist.begin(); it!= vlist.end(); it++ ) {
          cb[0] += it->E() ;
          cb[1] += it->Px() ;
          cb[2] += it->Py() ;
          cb[3] += it->Pz() ;
    }

    LorentzVector J12( cb[1], cb[2], cb[3], cb[0] );
    double mass2 =  (cb[0]*cb[0]) - (cb[1]*cb[1]) - (cb[2]*cb[2]) -  (cb[3]*cb[3]) ;
    if ( mass2 < 0 ) mass2 =0 ; 
    double mass = sqrt( mass2 );

    return mass;

}

double TtTools::getInvMass( LorentzVector lv ) {

     double mom2  = (lv.Px()*lv.Px()) + (lv.Py()*lv.Py()) + (lv.Pz()*lv.Pz()) ;
     double mass2 = lv.E()*lv.E() - mom2;
     double mass = 0. ;
     if (mass2 < 0. ) mass2 = 0;
     mass = sqrt(mass2) ;
     return mass;

}

double TtTools::getInvMass( LorentzVector lv1, LorentzVector lv2 ) {

    double p1 = lv1.Px() + lv2.Px() ;
    double p2 = lv1.Py() + lv2.Py() ;
    double p3 = lv1.Pz() + lv2.Pz() ;
    double p4 = lv1.E() + lv2.E() ;
    double mass2 = p4*p4 - (p1*p1) - (p2*p2) - (p3*p3) ;
    if (mass2 < 0. ) mass2 = 0;

    double  mass = sqrt(mass2) ;
    return mass;
}


double TtTools::getRelPt( LorentzVector a, LorentzVector b ) {

     // a x b = c
     double cx = (a.Py()*b.Pz()) - (a.Pz()*b.Py()) ;
     double cy = (a.Pz()*b.Px()) - (a.Px()*b.Pz()) ;
     double cz = (a.Px()*b.Py()) - (a.Py()*b.Px()) ;
     double cl = sqrt( cx*cx + cy*cy + cz*cz ) ;
     double ct = sqrt( cx*cx + cy*cy );
     double bl = sqrt( b.Px()*b.Px() + b.Py()*b.Py() + b.Pz()*b.Pz() );
     // a x b = |a|*|b|*sin(theta)
     double RelPt = (bl > 0 ) ? (cl/bl) : ct ;
     return RelPt ;

}

double TtTools::getBeta( LorentzVector a ) {

    double p = sqrt( a.Px()*a.Px() + a.Py()*a.Py() + a.Pz()*a.Pz() );
    double beta = ( a.E() <= 0.) ? -0.04 : p/a.E() ;
    if ( beta > 1. ) beta = 1.1 ;
    return beta ;

}

