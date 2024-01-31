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

using namespace o2;
using namespace o2::framework;

// quick task to read typical dileptonAod.root
// for now it only runs on the aod::DimuonsAll table but it can run on any table if you give it one

struct dileptonReader {

  // histograms can be created with OutputObj<TH1F>
  OutputObj<TH1F> hTauz{TH1F("Tauz", "Tauz", 100, -0.01, 0.01)};
  OutputObj<TH1F> hRapidityGlobal{TH1F("RapidityGlobal", "RapidityGlobal", 267, 2., 4.)};

  Configurable<Double_t> ptCut{"ptCut", 0.5, "pT cut, default = 0.5"}; // configurable pT cut

  // Create registry for histograms (several ways to create histograms)
  HistogramRegistry registry{
    "registry",
    {{"massGlobal", "M", {HistType::kTH1F, {{750, 0, 15}}}}}};

  void init(o2::framework::InitContext&)
  {
    // define axis for your histograms
    AxisSpec chi2Axis1 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 1"};
    AxisSpec chi2Axis2 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 2"};
    AxisSpec massAxis = {750, 0, 15, "M"};
    AxisSpec ptAxis = {120, 0.0, 30.0, "p_{T}"};

    HistogramConfigSpec ptSpec({HistType::kTH1F, {ptAxis}});
    HistogramConfigSpec histospec({HistType::kTH3F, {massAxis, chi2Axis1, chi2Axis2}});

    // add some histograms to the registry
    registry.add("Pt", "pt distribution", ptSpec);
    registry.add("massChi2TH3", "Mass and MCH-MFT Chi2 TH3 Histogram", histospec);
  };

  void process(aod::DimuonsAll const& dimuons)
  {
    // Double_t const PI = ROOT::Math::Pi(); // uncomment if you want to use pi, can be useful

    for (auto& dimuon : dimuons) {
      // calculate rapidity
      auto rap = -ROOT::Math::log((ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt() * ROOT::Math::cosh(dimuon.eta()) * ROOT::Math::cosh(dimuon.eta())) + dimuon.pt() * ROOT::Math::sinh(dimuon.eta())) / (ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt())));

      // MUON standalone tracks studies
      if ((dimuon.eta1() > -4 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -4 && dimuon.eta2() < -2.5)) {
        if (dimuon.pt1() >= ptCut && dimuon.pt2() >= ptCut) {
          if (dimuon.sign() == 0) {
            if (dimuon.chi2MatchMCHMFT1() <= 0. && dimuon.chi2MatchMCHMFT2() <= 0.) {
              if (rap > 2.5 && rap < 4) {
                // fill any histo for MCH-MID tracks (create them first in the histogram registry OR define them befor Init() with OutputObj)
              } // end rapidity cut
            }   // end MCH-MID track selection
          }     // end OS selection
        }       // end pt cut
      }         // end eta cut (MUON standalone acceptance)

      // Global tracks studies
      if ((dimuon.eta1() > -3.6 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -3.6 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
        if (dimuon.chi2MatchMCHMFT1() > 0 && dimuon.chi2MatchMCHMFT2() > 0) {
          if (!(dimuon.isAmbig1()) && !(dimuon.isAmbig2())) { // remove ambiguous tracks
            hRapidityGlobal->Fill(rap);                       // fill rapidity distribution histogram
            if (dimuon.sign() == 0) {
              if (rap > 2.5 && rap < 3.6) {
                hTauz->Fill(dimuon.tauz());

                registry.get<TH1>(HIST("massGlobal"))->Fill(dimuon.mass());
                registry.get<TH1>(HIST("Pt"))->Fill(dimuon.pt());

                // TH3 containing mass and chi2 information filled to do chi2 study                                                    // fill mass invariant distribution for Global tracks
                registry.get<TH3>(HIST("massChi2TH3"))->Fill(dimuon.mass(), dimuon.chi2MatchMCHMFT1(), dimuon.chi2MatchMCHMFT2());

              } // end rapidity cut
            }   // end sign selection
          }     // ambig cut
        }       // end MCH-MFT chi2 selection
      }         // end eta cut
    }           // end loop over dimuons
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dileptonReader>(cfgc)};
};
