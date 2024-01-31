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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"

#include <TFile.h>
#include <TObject.h>
#include <TDirectoryFile.h>
#include <TH1F.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;


using MyMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

struct muonLoopEventReader {

  // histogram registry delclaration
  HistogramRegistry registry{
    "registry",
    {{"mass", "M", {HistType::kTH1F, {{750, 0, 15}}}}}};

  void init(o2::framework::InitContext&)
  {
    // Event histograms
    AxisSpec multAxis = {100, 0., 100.};
    HistogramConfigSpec multSpec({HistType::kTH1F, {multAxis}});
    registry.add("multNContribPV", "multiplcity (N contributor to PV)", multSpec);


    AxisSpec massAxis = {750, 0, 15, "m_{#mu#mu}"};
    AxisSpec ptAxis = {200, 0.0, 20.0, "p_{T}"};
    AxisSpec pxAxis = {200, -10.0, 10.0, "p_{x}"};
    AxisSpec pyAxis = {200, -10.0, 10.0, "p_{y}"};
    AxisSpec rapAxis = {200, 2.5, 4.0, "y"};
    AxisSpec chi2Axis = {250, 0., 500., "#chi^{2}_{MCH-MFT}"};
    AxisSpec etaAxis = {500, -5.0, 5.0, "#eta"};
    AxisSpec rAbsAxis = {1000, 0., 100., "R_{abs}"};

    HistogramConfigSpec ptSpec({HistType::kTH1F, {ptAxis}});
    HistogramConfigSpec pxSpec({HistType::kTH1F, {pxAxis}});
    HistogramConfigSpec pySpec({HistType::kTH1F, {pyAxis}});
    HistogramConfigSpec etaSpec({HistType::kTH1F, {etaAxis}});
    HistogramConfigSpec chi2Spec({HistType::kTH1F, {chi2Axis}});
    HistogramConfigSpec rAbsSpec({HistType::kTH1F, {rAbsAxis}});
    HistogramConfigSpec chi2vsPtSpec({HistType::kTH2F, {ptAxis, chi2Axis}});
    HistogramConfigSpec chi2VsPtVsRabsSpec({HistType::kTH3F, {ptAxis, chi2Axis, rAbsAxis}});
    HistogramConfigSpec chi2VsPtVsNContribPVSpec({HistType::kTH3F, {ptAxis, chi2Axis, multAxis}});
    HistogramConfigSpec massSpec({HistType::kTH1F, {massAxis}});

    registry.add("chi2MCHMFT", "chi2MCHMFT", chi2Spec);
    registry.add("chi2MCHMFT_vs_pT", "Chi2MCHMFT vs pT", chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_rAbs_23", "Chi2MCHMFT vs pT with 2 < #theta < 3", chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_rAbs_310", "Chi2MCHMFT vs pT with 3 < #theta < 10", chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_vs_rAbs", "Chi2MCHMFT vs pT vs R_{abs}", chi2VsPtVsRabsSpec);
    registry.add("chi2MCHMFT_vs_pT_left", "chi2MCHMFT vs pT left",chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_right", "chi2MCHMFT vs pT right",chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_top", "chi2MCHMFT vs pT top",chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_bottom", "chi2MCHMFT vs pT bottom",chi2vsPtSpec);
    registry.add("chi2MCHMFT_vs_pT_vs_NContribPV", "chi2MCHMFT vs pT vs N contrib PV",chi2VsPtVsNContribPVSpec);

    registry.add("pt", "p_{#rm{T}}", ptSpec);
    registry.add("px", "p_{x}", pxSpec);
    registry.add("py", "p_{y}", pySpec);
    registry.add("eta", "#eta", etaSpec);
    registry.add("rAbs", "R_{abs}", rAbsSpec);

    registry.add("chi2MCHMFT_withAmbiguous", "chi2MCHMFT_withAmbiguous", chi2Spec);
    registry.add("chi2MCHMFT_vs_pT_withAmbiguous", "Chi2MCHMFT vs pT_withAmbiguous", chi2vsPtSpec);
    registry.add("pt_withAmbiguous", "p_{#rm{T}}_withAmbiguous", ptSpec);
    registry.add("eta_withAmbiguous", "#eta_withAmbiguous", etaSpec);

    // temp test
    registry.add("Mass", "m_{#mu#mu}", massSpec);
    
  }

  void process(MyEvents::iterator const& event, MyMuons const& muons)
  {
    registry.get<TH1>(HIST("multNContribPV"))->Fill(event.multNTracksPV());
    for (auto& muon : muons) {    
      registry.get<TH1>(HIST("pt"))->Fill(muon.pt());
      registry.get<TH3>(HIST("chi2MCHMFT_vs_pT_vs_NContribPV"))->Fill(muon.pt(), muon.chi2MatchMCHMFT(), event.multNTracksPV());
    }
    
  }
};

// struct muonReaderMC {
//   OutputObj<TH1F> hChi2MCHMFT{TH1F("Chi2MCHMFT", "Chi2MCHMFT", 100, 0., 200.)};
//   OutputObj<TH1F> hChi2MCHMFTTrue{TH1F("Chi2MCHMFTTrue", "Chi2MCHMFTTrue", 100, 0., 200.)};
//   OutputObj<TH1F> hChi2MCHMFTFalse{TH1F("Chi2MCHMFTFalse", "Chi2MCHMFTFalse", 100, 0., 200.)};

//   void init(o2::framework::InitContext&)
//   {
//   }

//   void process(MyMuonsMC const& muons)
//   {
//     for (auto& muon : muons) {
//       hChi2MCHMFT->Fill(muon.chi2MatchMCHMFT());
//       if (muon.mcMask() < 20) {
//         hChi2MCHMFTTrue->Fill(muon.chi2MatchMCHMFT());
//       }
//       if (muon.mcMask() > 20) {
//         hChi2MCHMFTFalse->Fill(muon.chi2MatchMCHMFT());
//       }
//     }
//   }
// };

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    // adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<muonLoopEventReader>(cfgc)
    // adaptAnalysisTask<muonReaderMC>(cfgc),
  };
};