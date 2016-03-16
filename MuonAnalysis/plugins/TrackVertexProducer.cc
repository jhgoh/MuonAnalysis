#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

using namespace std;

class TrackVertexProducer : public edm::stream::EDProducer<>
{
public:
  TrackVertexProducer(const edm::ParameterSet& pset);
  virtual ~TrackVertexProducer() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  std::pair<reco::VertexCompositeCandidate, double> fitSV(const reco::Vertex& pv,
                                                          const reco::TransientTrack& transTrack1,
                                                          const reco::TransientTrack& transTrack2) const;

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  const double pionMass = 0.1396;
  const double kaonMass = 0.4937;
  const double protonMass = 0.9383;

  const double trkMinPt_, trkMaxEta_;
  double vtxMinRawMass_, vtxMaxRawMass_;
  double vtxMinMass_, vtxMaxMass_, vtxMinLxy_, vtxMaxLxy_;
  const double trkChi2_, trkDCA_, vtxChi2_, vtxSignif_;

  int vtxPdgId_;
  double mass1_, mass2_;
  bool isSameFlav_;
};

TrackVertexProducer::TrackVertexProducer(const edm::ParameterSet& pset):
  trkMinPt_(pset.getParameter<double>("trkMinPt")),
  trkMaxEta_(pset.getParameter<double>("trkMaxEta")),
  vtxMinLxy_(pset.getParameter<double>("vtxMinLxy")),
  vtxMaxLxy_(pset.getParameter<double>("vtxMaxLxy")),
  trkChi2_(pset.getParameter<double>("trkChi2")),
  trkDCA_(pset.getParameter<double>("trkDCA")),
  vtxChi2_(pset.getParameter<double>("vtxChi2")),
  vtxSignif_(pset.getParameter<double>("vtxSignif"))
{
  trackToken_ = consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("src"));
  vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));

  const string vtxType = pset.getParameter<string>("vtxType");
  if ( vtxType == "kshort" ) {
    vtxPdgId_ = 310;
    isSameFlav_ = true;
    mass1_ = pionMass; mass2_ = pionMass;
    vtxMinRawMass_ = 0.40; vtxMaxRawMass_ = 0.60;
    vtxMinMass_ = 0.43; vtxMaxMass_ = 0.57;
  }
  else if ( vtxType == "phi" ) {
    vtxPdgId_ = 333;
    isSameFlav_ = true;
    mass1_ = kaonMass; mass2_ = kaonMass;
    vtxMinRawMass_ = 0.95; vtxMaxRawMass_ = 1.08;
    vtxMinMass_ = 0.98; vtxMaxMass_ = 1.06;
  }
  else if ( vtxType == "lambda" ) {
    vtxPdgId_ = 3122;
    mass1_ = protonMass; mass2_ = pionMass;
    vtxMinRawMass_ = 1.091; vtxMaxRawMass_ = 1.139;
    vtxMinMass_ = 1.101; vtxMaxMass_ = 1.129;
  }
  else {
    vtxPdgId_ = 0;
    mass1_ = pset.getParameter<double>("mass1");
    mass2_ = pset.getParameter<double>("mass2");
    isSameFlav_ = (mass1_ == mass2_);
    vtxMinRawMass_ = pset.getParameter<double>("minRawMass");
    vtxMaxRawMass_ = pset.getParameter<double>("maxRawMass");
    vtxMinMass_ = pset.getParameter<double>("minMass");
    vtxMaxMass_ = pset.getParameter<double>("maxMass");
  }

  produces<reco::VertexCompositeCandidateCollection>();
  produces<std::vector<double> >("lxy");
}

void TrackVertexProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<reco::VertexCompositeCandidateCollection> out(new reco::VertexCompositeCandidateCollection);
  std::auto_ptr<std::vector<double> > out_lxy(new std::vector<double>());

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  edm::Handle<reco::TrackCollection> trackHandle;
  event.getByToken(trackToken_, trackHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  reco::Vertex pv = vertexHandle->at(0);

  std::vector<reco::TransientTrack> transTracks;
  for ( auto track = trackHandle->begin(); track != trackHandle->end(); ++track ) {
    if ( track->pt() < trkMinPt_ or std::abs(track->eta()) > trkMaxEta_ ) continue;
    // Apply basic track quality cuts
    auto transTrack = trackBuilder->build(&*track);
    transTracks.push_back(transTrack);
  }

  for ( auto itr1 = transTracks.begin(); itr1 != transTracks.end(); ++itr1 ) {
    const reco::Track& track1 = itr1->track();
    const double e1 = sqrt(mass1_*mass1_ + track1.momentum().mag2());

    for ( auto itr2 = isSameFlav_ ? itr1+1 : transTracks.begin(); itr2 != transTracks.end(); ++itr2 ) {
      if ( itr1 == itr2 ) continue;
      if ( itr1->charge() == itr2->charge() ) continue;
      const reco::Track& track2 = itr2->track();
      const double e2 = sqrt(mass2_*mass2_ + track2.momentum().mag2());

      const double px = track1.px() + track2.px();
      const double py = track1.px() + track2.px();
      const double pz = track1.px() + track2.px();
      const double p2 = px*px + py*py + pz*pz;
      const double e = e1+e2;

      const double rawMass = sqrt(e*e - p2);
      if ( rawMass < vtxMinRawMass_ or rawMass > vtxMaxRawMass_ ) continue;

      auto sv = fitSV(pv, *itr1, *itr2);
      if ( sv.first.pdgId() == 0 ) continue;

      out->push_back(sv.first);
      out_lxy->push_back(sv.second);
    }
  }

  event.put(out);
  event.put(out_lxy, "lxy");
}

std::pair<reco::VertexCompositeCandidate, double> TrackVertexProducer::fitSV(const reco::Vertex& pv,
                                                                             const reco::TransientTrack& transTrack1,
                                                                             const reco::TransientTrack& transTrack2) const
{
  std::pair<reco::VertexCompositeCandidate, double> result;
  result.first.setPdgId(0);
  result.second = -1;

  try {
    auto ipState1 = transTrack1.impactPointTSCP().theState();
    auto ipState2 = transTrack2.impactPointTSCP().theState();

    ClosestApproachInRPhi cApp;
    cApp.calculate(ipState1, ipState2);
    if ( !cApp.status() ) return result;

    const float dca = std::abs(cApp.distance());
    if ( dca < 0. or dca > trkDCA_ ) return result;

    GlobalPoint cxPt = cApp.crossingPoint();
    if ( std::hypot(cxPt.x(), cxPt.y()) > 120. or std::abs(cxPt.z()) > 300. ) return result;

    //TrajectoryStateClosestToPoint caState1 = leptonTRack.trajectoryStateClosestToPoint(cxPt);
    //TrajectoryStateClosestToPoint caState2 = transTrackNeg.trajectoryStateClosestToPoint(cxPt);
    //if ( !caState1.isValid() or !caState2.isValid() ) return result;

    std::vector<reco::TransientTrack> transTracks = {transTrack1, transTrack2};

    KalmanVertexFitter fitter(true);
    TransientVertex tsv = fitter.vertex(transTracks);
    if ( !tsv.isValid() or tsv.totalChiSquared() < 0. ) return result;

    reco::Vertex sv = tsv;
    if ( sv.normalizedChi2() > vtxChi2_ ) return result;

    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
    typedef ROOT::Math::SVector<double, 3> SVector3;

    GlobalPoint vtxPos(sv.x(), sv.y(), sv.z());
    SMatrixSym3D totalCov = pv.covariance() + sv.covariance();
    SVector3 distanceVectorXY(sv.x() - pv.x(), sv.y() - pv.y(), 0.);

    const double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
    const double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;
    if( rVtxMag < vtxMinLxy_ or rVtxMag > vtxMaxLxy_ or rVtxMag / sigmaRvtxMag < vtxSignif_ ) return result;

    //SVector3 distanceVector3D(sv.x() - pvx, sv.y() - pvy, sv.z() - pvz);
    //const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);

    // Cuts finished, now we create the candidates and push them back into the collections.
    GlobalVector mom1, mom2;
    if ( !tsv.hasRefittedTracks() ) {
      mom1 = transTrack1.trajectoryStateClosestToPoint(vtxPos).momentum();
      mom2 = transTrack2.trajectoryStateClosestToPoint(vtxPos).momentum();
    }
    else {
      auto refTracks = tsv.refittedTracks();
      if ( refTracks.size() < 2 ) return result;

      mom1 = refTracks.at(0).trajectoryStateClosestToPoint(vtxPos).momentum();
      mom2 = refTracks.at(1).trajectoryStateClosestToPoint(vtxPos).momentum();
    }
    if ( mom1.mag() <= 0 or mom2.mag() <= 0 ) return result;
    const GlobalVector mom = mom1+mom2;

    const double candE1 = hypot(mom1.mag(), mass1_);
    const double candE2 = hypot(mom2.mag(), mass2_);
    const double vtxChi2 = sv.chi2();
    const double vtxNdof = sv.ndof();

    reco::Particle::Point vtx(sv.x(), sv.y(), sv.z());
    const reco::Vertex::CovarianceMatrix vtxCov(sv.covariance());

    const math::XYZTLorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
    if ( vtxMinMass_ > candLVec.mass() or vtxMaxMass_ < candLVec.mass() ) return result;

    result.first.setPdgId(vtxPdgId_);
    result.first.setP4(candLVec);
    result.first.setVertex(vtx);
    result.first.setChi2AndNdof(vtxChi2, vtxNdof);
    result.first.setCovariance(vtxCov);
    result.second = rVtxMag;

  } catch ( std::exception& e ) { return result; }

  return result;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackVertexProducer);

