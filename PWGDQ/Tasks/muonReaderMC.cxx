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
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"

#include <TFile.h>
#include <TObject.h>
#include <TDirectoryFile.h>
#include <TH1F.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;

using MyMuonsMC = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;

struct muonReaderMC {
  // histogram registry delclaration
  HistogramRegistry registry{
    "registry",
    {{"mass", "M", {HistType::kTH1F, {{750, 0, 15}}}}}};

  void init(o2::framework::InitContext&)
  {
    AxisSpec massAxis = {750, 0, 15, "M"};
    AxisSpec ptAxis = {2000, 0.0, 20.0, "p_{T}"};
    AxisSpec rapAxis = {200, 2.5, 4.0, "y"};
    AxisSpec chi2Axis = {250, 0., 500., "#chi^{2}_{MCH-MFT}"};
    AxisSpec etaAxis = {500, -5.0, 5.0, "#eta"};
    AxisSpec maskAxis = {200, 0, 200, "mcMask"};
    

    HistogramConfigSpec ptSpec({HistType::kTH1F,{ptAxis}});
    HistogramConfigSpec etaSpec({HistType::kTH1F,{etaAxis}});
    HistogramConfigSpec chi2Spec({HistType::kTH1F,{chi2Axis}});
    HistogramConfigSpec chi2vsPtSpec({HistType::kTH2F, {ptAxis, chi2Axis}});
    HistogramConfigSpec chi2vsPtvsMcMaskSpec({HistType::kTH3F, {chi2Axis,ptAxis,maskAxis}});

    registry.add("MCchi2MCHMFT", "chi2MCHMFT", chi2Spec);
    registry.add("MCchi2MCHMFTTrue", "chi2MCHMFT true matches", chi2Spec);
    registry.add("MCchi2MCHMFTFalse", "chi2MCHMFT false matches", chi2Spec);

    registry.add("MCchi2MCHMFT_vs_pT", "Chi2MCHMFT vs pT", chi2vsPtSpec);
    registry.add("MCchi2MCHMFT_vs_pT_true", "Chi2MCHMFT vs pT - true matches", chi2vsPtSpec);
    registry.add("MCchi2MCHMFT_vs_pT_false", "Chi2MCHMFT vs pT - false matches", chi2vsPtSpec);

    registry.add("MCTH3Chi2PtMcMask", "TH3 chi2 vs pT vs McMask", chi2vsPtvsMcMaskSpec);
    registry.add("MCpt","p_{#rm{T}}",ptSpec);
    registry.add("MCeta","#eta",etaSpec);
  }

  void process(MyMuonsMC const& muons)
  {
    for (auto& muon : muons) {

      registry.get<TH1>(HIST("MCpt"))->Fill(muon.pt());
      registry.get<TH1>(HIST("MCeta"))->Fill(muon.eta());
      registry.get<TH1>(HIST("MCchi2MCHMFT"))->Fill(muon.chi2MatchMCHMFT());
      registry.get<TH2>(HIST("MCchi2MCHMFT_vs_pT"))->Fill(muon.pt(), muon.chi2MatchMCHMFT());
      
      registry.get<TH3>(HIST("MCTH3Chi2PtMcMask"))->Fill(muon.chi2MatchMCHMFT(),muon.pt(),muon.mcMask());
      
      if (muon.mcMask() < 20) {
        registry.get<TH1>(HIST("MCchi2MCHMFTTrue"))->Fill(muon.chi2MatchMCHMFT());
        registry.get<TH2>(HIST("MCchi2MCHMFT_vs_pT_true"))->Fill(muon.pt(), muon.chi2MatchMCHMFT());
      }
      if (muon.mcMask() > 20) {
        registry.get<TH1>(HIST("MCchi2MCHMFTFalse"))->Fill(muon.chi2MatchMCHMFT());
        registry.get<TH2>(HIST("MCchi2MCHMFT_vs_pT_false"))->Fill(muon.pt(), muon.chi2MatchMCHMFT());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<muonReaderMC>(cfgc)
    };
};
