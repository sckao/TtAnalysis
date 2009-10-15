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
  double al = sqrt( (v1.Px()*v1.Px()) + (v1.Py()*v1.Py()) );
  double bl = sqrt( (v2.Px()*v2.Px()) + (v2.Py()*v2.Py()) );
  double cosA = ab/(al*bl) ;

  double dv = ab - al*bl ;
  if ( dv > 0. && dv < 0.0001 ) cosA = 1. ;
  double df = acos(cosA) ;

  bool normal = ( df >= 0. && df < 3.1415927 ) ? true : false ;
  if ( !normal ) cout<<" df:"<<df<<" v1:"<<al<<" v2:"<<bl<<"  v1v2:"<<ab<<endl;

  return df;
}

double TtTools::get_dPhi( LorentzVector v1, LorentzVector v2 ) {

  double axb = (v1.Px()*v2.Py()) - (v1.Py()*v2.Px()) ;
  double ab = (v1.Px()*v2.Px()) + (v1.Py()*v2.Py()) ;
  double al = sqrt( (v1.Px()*v1.Px()) + (v1.Py()*v1.Py()) );
  double bl = sqrt( (v2.Px()*v2.Px()) + (v2.Py()*v2.Py()) );
  double cosA = ab/(al*bl) ;

  double df = acos(cosA) ;
  if ( axb < 0 ) df = -1.*df ;

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

  //if ( dR < 0.000001 ) dR = 0.0 ;
  bool normal = ( dR <  9. ) ? true:false ;
  if ( !normal ) cout <<" Y1:"<< aY <<" Y2:"<<bY<<" df:"<<df<< endl;

  return dR;
}

double TtTools::getY( LorentzVector v1 ){

    double ep = (v1.E() + v1.Pz()) / (v1.E() - v1.Pz()) ;
    double Y = 0.5*log( ep  ) ;
    if ( (v1.E() - v1.Pz()) == 0. ) Y =  99.99;
    if ( (v1.E() + v1.Pz()) == 0. ) Y = -99.99;
    return Y;
}

// for leptonic W
double TtTools::getMt( LorentzVector v1,  LorentzVector v2 ) {

     double Et = ( v1.Et() + v2.Et() ) ;
     double px = ( v1.Px() + v2.Px() ) ;
     double py = ( v1.Py() + v2.Py() ) ;
     double Pt = sqrt( (px*px) + (py*py) ) ;
     double Mt2 = (Et*Et) - (Pt*Pt) ;
     double Mt = ( Mt2 > 0. ) ? sqrt(Mt2) : 0 ;

     return Mt ;
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

std::vector<const pat::Jet*> TtTools::ReturnJetForm( std::vector<const reco::Candidate*> jCand, Handle<std::vector<pat::Jet> > patJet ){

     std::vector<const pat::Jet*> emptyCont ;
     std::vector<const pat::Jet*> jet_temp ;
     std::vector<int> idx ;
     bool fail = false ;
     for (size_t i=0; i< jCand.size(); i++ ) {
         int jdx = -1;
         for ( size_t j=0; j<  patJet->size(); j++) {
             double dx = jCand[i]->px() - (*patJet)[j].px() ;
             double dy = jCand[i]->py() - (*patJet)[j].py() ;
             double dz = jCand[i]->pz() - (*patJet)[j].pz() ;
             double de = jCand[i]->energy() - (*patJet)[j].energy() ;
             double dAll = sqrt( dx*dx + dy*dy + dz*dz + de*de ) ;
             if ( dAll < 0.00001 ) {
                jet_temp.push_back( &(*patJet)[j] );
                jdx = static_cast<int>(j);
             }
         }
         for (size_t k=0; k< idx.size(); k++ ) {
             if (idx[k] == jdx ) fail = true ;
         }
         if ( jdx != -1 && !fail )  idx.push_back( jdx );
     }
     if ( jet_temp.size() != jCand.size() ) return emptyCont;
     if ( fail ) return emptyCont;
     return jet_temp;
}

// match the reco candidate with pat jets
const pat::Jet* TtTools::ReturnJetForm( const reco::Candidate* jCand, Handle<std::vector<pat::Jet> > patJet, bool& goodmatching ){

     const pat::Jet* patform = NULL;
     int idx = -1 ;
     double dAll0 = 0.1;
     for ( size_t j=0; j<  patJet->size(); j++) {
         double dx = jCand->px() - (*patJet)[j].px() ;
         double dy = jCand->py() - (*patJet)[j].py() ;
         double dz = jCand->pz() - (*patJet)[j].pz() ;
         double de = jCand->energy() - (*patJet)[j].energy() ;
         double dAll = sqrt( dx*dx + dy*dy + dz*dz + de*de ) ;
         if ( dAll < dAll0 ) {
            idx = static_cast<int>(j);
            dAll0 = dAll;
         }
     }
     if ( idx != -1 ) patform =  &(*patJet)[idx] ;
     if ( dAll0 < 0.0001 ) goodmatching = true; 

     return patform ;
}
