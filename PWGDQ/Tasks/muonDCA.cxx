// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \file muonDCA.cxx
/// \brief Task to compute and evaluate DCA quantities
/// \author Nicolas Bizé <nicolas.bize@cern.ch>, SUBATECH

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

static o2::globaltracking::MatchGlobalFwd mExtrap;
template <typename T>
bool isSelected(const T& muon);
template <typename T, typename C>
o2::dataformats::GlobalFwdTrack propagateToVertex(const T& muon, const C& collision);
template <typename T, typename C>
o2::dataformats::GlobalFwdTrack propagateToDCA(const T& muon, const C& collision);
template <typename T, typename C>
o2::dataformats::GlobalFwdTrack propagateToRabs(const T& muon, const C& collision);

struct muonExtrap {
  Produces<ReducedMuonsDca> dcaTable;
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::parameters::GRPMagField* grpmag = nullptr; // for run 3, we access GRPMagField from GLO/Config/GRPMagField
  int fCurrentRun;                               // needed to detect if the run changed and trigger update of magnetic field

  HistogramRegistry registry{
    "registry",
    {}};

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      fCCDB->get<TGeoManager>(geoPath);
    }

    AxisSpec pdcaAxis = {5000, 0.0, 5000.0, "p #times DCA"};
    AxisSpec dcaAxis = {200, 0.0, 200.0, "DCA"};
    AxisSpec dcaxAxis = {200, -100.0, 100.0, "DCA_x"};
    AxisSpec dcayAxis = {200, -100.0, 100.0, "DCA_y"};
    AxisSpec rabsAxis = {100, 0., 100.0, "R_{abs}"};

    HistogramConfigSpec pdcaSpec({HistType::kTH1F, {pdcaAxis}});
    HistogramConfigSpec dcaSpec({HistType::kTH1F, {dcaAxis}});
    HistogramConfigSpec dcaxSpec({HistType::kTH1F, {dcaxAxis}});
    HistogramConfigSpec dcaySpec({HistType::kTH1F, {dcayAxis}});
    HistogramConfigSpec rabsSpec({HistType::kTH1F, {rabsAxis}});

    registry.add("pdca", "pDCA", pdcaSpec);
    registry.add("dca", "DCA", dcaSpec);
    registry.add("dcax", "DCA_x", dcaxSpec);
    registry.add("dcay", "DCA_y", dcaySpec);
    registry.add("rabs", "R_{abs}", rabsSpec);
  }

  void processExtrapolation(MyEventsVtxCov::iterator const& collision, MyMuons const& muons)
  {
    if (fCurrentRun != collision.runNumber()) {
      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, collision.timestamp());
      if (grpmag != nullptr) {
        LOGF(info, "Init field from GRP");
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
      LOGF(info, "Set field for muons");
      VarManager::SetupMuonMagField();
      fCurrentRun = collision.runNumber();
    }

    for (auto& muon : muons) {
      if (static_cast<int>(muon.trackType()) < 2) {
        continue; // Make sure to remove global muon tracks
      }
      // propagate muon track to vertex
      o2::dataformats::GlobalFwdTrack muonTrackAtVertex = propagateToVertex(muon, collision);

      // propagate muon track to DCA
      o2::dataformats::GlobalFwdTrack muonTrackAtDCA = propagateToDCA(muon, collision);

      // propagate to Rabs
      o2::dataformats::GlobalFwdTrack muonTrackAtRabs = propagateToRabs(muon, collision);

      // Calculate DCA quantities (preferable to do it with VarManager)
      double dcax = collision.posX() - muonTrackAtDCA.getX();
      double dcay = collision.posY() - muonTrackAtDCA.getY();
      double dca = std::sqrt(dcax * dcax + dcay * dcay);
      double pdca = muonTrackAtVertex.getP() * dca;

      double xAbs = muonTrackAtRabs.getX();
      double yAbs = muonTrackAtRabs.getY();
      double rabs = std::sqrt(xAbs * xAbs + yAbs * yAbs);

      // QA histograms
      registry.get<TH1>(HIST("pdca"))->Fill(pdca);
      registry.get<TH1>(HIST("dca"))->Fill(dca);
      registry.get<TH1>(HIST("dcax"))->Fill(dcax);
      registry.get<TH1>(HIST("dcay"))->Fill(dcay);
      registry.get<TH1>(HIST("rabs"))->Fill(rabs);

      // Fill DCA table
      dcaTable(pdca,
               dca,
               dcax,
               dcay,
               rabs,
               muonTrackAtVertex.getPt(),
               muonTrackAtVertex.getEta(),
               muonTrackAtVertex.getPhi(),
               muon.sign(),
               muon.isAmbiguous());
    }
  }

  PROCESS_SWITCH(muonExtrap, processExtrapolation, "process extrapolation", false);

  void processDummy(MyEventsVtxCov&)
  {
    // do nothing
  }

  PROCESS_SWITCH(muonExtrap, processDummy, "do nothing", false);
};

template <typename T>
bool isSelected(const T& muon)
{
  bool keepTrack = true;
  if (muon.eta() < -4. || muon.eta() > -2.5) {
    LOGF(info, "Reject muon with eta = %f", muon.eta());
    keepTrack = false;
  }
  if (muon.isAmbiguous()) {
    LOGF(info, "Reject ambiguous muon with flag = %d", muon.isAmbiguous());
    keepTrack = false;
  }
  return keepTrack;
}


// template <typename T, typename C>
// o2::dataformats::GlobalFwdTrack propagateGeneralMethod(const T& muon, const C& collision)
// {
//   double chi2 = muon.chi2();
//   SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
//   std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
//                          muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
//                          muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
//   SMatrix55 tcovs(v1.begin(), v1.end());
//   o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
//   o2::dataformats::GlobalFwdTrack propmuon;
//   if (static_cast<int>(muon.trackType()) > 2) {
//     o2::dataformats::GlobalFwdTrack track;
//     track.setParameters(tpars);
//     track.setZ(fwdtrack.getZ());
//     track.setCovariances(tcovs);
//     auto mchTrack = mExtrap.FwdtoMCH(track);

//     if (kToVertex){
//       // do vertex
//       o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
//     auto proptrack = mExtrap.MCHtoFwd(mchTrack);
//     propmuon.setParameters(proptrack.getParameters());
//     propmuon.setZ(proptrack.getZ());
//     propmuon.setCovariances(proptrack.getCovariances());
//     return propmuon;
//     }
//     else if (kToDCA){
//       // do dca
//     o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
//     auto proptrack = mExtrap.MCHtoFwd(mchTrack);
//     propmuon.setParameters(proptrack.getParameters());
//     propmuon.setZ(proptrack.getZ());
//     propmuon.setCovariances(proptrack.getCovariances());
//     return propmuon;
//     }
//     else if (kToRabs){
//       // do rabs
//     o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
//     auto proptrack = mExtrap.MCHtoFwd(mchTrack);
//     propmuon.setParameters(proptrack.getParameters());
//     propmuon.setZ(proptrack.getZ());
//     propmuon.setCovariances(proptrack.getCovariances());
//     return propmuon;
//     }
//     o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
//     auto proptrack = mExtrap.MCHtoFwd(mchTrack);
//     propmuon.setParameters(proptrack.getParameters());
//     propmuon.setZ(proptrack.getZ());
//     propmuon.setCovariances(proptrack.getCovariances());
//   }
//   return propmuon;
// }

template <typename T, typename C>
o2::dataformats::GlobalFwdTrack propagateToVertex(const T& muon, const C& collision)
{
  double chi2 = muon.chi2();
  SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
  std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                         muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                         muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
  o2::dataformats::GlobalFwdTrack propmuon;
  if (static_cast<int>(muon.trackType()) > 2) {
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = mExtrap.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    auto proptrack = mExtrap.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
  }
  return propmuon;
}

template <typename T, typename C>
o2::dataformats::GlobalFwdTrack propagateToDCA(const T& muon, const C& collision)
{
  double chi2 = muon.chi2();
  SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
  std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                         muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                         muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
  o2::dataformats::GlobalFwdTrack propmuon;
  if (static_cast<int>(muon.trackType()) > 2) {
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = mExtrap.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
    auto proptrack = mExtrap.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
  }
  return propmuon;
}

template <typename T, typename C>
o2::dataformats::GlobalFwdTrack propagateToRabs(const T& muon, const C& collision)
{
  double chi2 = muon.chi2();
  SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
  std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                         muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                         muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
  o2::dataformats::GlobalFwdTrack propmuon;
  // o2::mch::TrackParam mchTrack;
  if (static_cast<int>(muon.trackType()) > 2) {
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = mExtrap.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
    auto proptrack = mExtrap.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
  }
  return propmuon;
}

// extrapolate to the end of the absorber
// o2::mch::TrackParam trackParamAtRAbs(track.getZ(), track.getParameters());
// if (!o2::mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
//   return false;
// }

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<muonExtrap>(cfgc)};
};