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

struct mumuReader {

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> hchi2MCHMIDLeft{TH1F("chi2MCHMIDLeft", "chi2MCHMIDLeft", 100, 0., 16.)};
  OutputObj<TH1F> hchi2MCHMIDRight{TH1F("chi2MCHMIDRight", "chi2MCHMIDRight", 100, 0., 16.)};
  OutputObj<TH1F> hchi2MCHMIDUp{TH1F("chi2MCHMIDUp", "chi2MCHMIDUp", 100, 0., 16.)};
  OutputObj<TH1F> hchi2MCHMIDDown{TH1F("chi2MCHMIDDown", "chi2MCHMIDDown", 100, 0., 16.)};
  // check MCH-MFT chi2 also
  OutputObj<TH1F> hChi2MCHMFTLeft{TH1F("chi2MCHMFTLeft", "chi2MCHMFTLeft", 100, 0., 16.)};
  OutputObj<TH1F> hChi2MCHMFTRight{TH1F("chi2MCHMFTRight", "chi2MCHMFTRight", 100, 0., 16.)};

  OutputObj<TH1F> hRapidityAll{TH1F("Rapidity", "Rapidity", 267, 2., 4.)};

  OutputObj<TH1F> hRapidityMCHMFT{TH1F("RapidityMCHMFT", "RapidityMCHMFT", 267, 2., 4.)};
  OutputObj<TH1F> hTauz{TH1F("Tauz", "Tauz", 100, -0.01, 0.01)};
  OutputObj<TH1F> hdeltaZ{TH1F("deltaZ", "deltaZ", 1000, -10., 10.)};
  OutputObj<TH1F> hChi2MCHMFT{TH1F("Chi2MCHMFT", "Chi2MCHMFT", 250, -50., 200.)};

  Configurable<bool> usePhi{"usePhi", false, "If true, use phi method"};
  Configurable<bool> useMomentum{"useMomentum", false, "If true, use momentum method"};
  Configurable<Double_t> ptCut{"ptCut", 0.5, "pt cut, default = 0.5"};

  // histogram registry test
  HistogramRegistry registry{
    "registry",
    {{"mass", "M", {HistType::kTH1F, {{750, 0, 15}}}}}};
  // end of histogram registry test

  void init(o2::framework::InitContext&)
  {
    AxisSpec massAxis = {750, 0, 15, "M"};
    AxisSpec chi2Axis1 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 1"};
    AxisSpec chi2Axis2 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 2"};
    HistogramConfigSpec histospec({HistType::kTH3F, {massAxis, chi2Axis1, chi2Axis2}});
    HistogramConfigSpec massSpec({HistType::kTH1F, {massAxis}});

    // MUON standalone studies
    registry.add("JPsiXPos", "JPsi X pos", massSpec);
    registry.add("JPsiXNeg", "JPsi X neg", massSpec);
    registry.add("JPsiXNeutral", "JPsi X neutral", massSpec);
    registry.add("JPsiYPos", "JPsi Y pos", massSpec);
    registry.add("JPsiYNeg", "JPsi Y neg", massSpec);
    registry.add("JPsiYNeutral", "JPsi Y neutral", massSpec);
    registry.add("JPsiStandard", "JPsi standard", massSpec);

    // MCH-MFT chi2 studies
    registry.add("massTH3", "Mass TH3 Histogram", histospec);
  };

  void process(aod::DimuonsAll const& dimuons)
  {
    Double_t const PI = ROOT::Math::Pi();

    for (auto& dimuon : dimuons) {
      // dimuon variables
      auto rap = -ROOT::Math::log((ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt() * ROOT::Math::cosh(dimuon.eta()) * ROOT::Math::cosh(dimuon.eta())) + dimuon.pt() * ROOT::Math::sinh(dimuon.eta())) / (ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt())));
      auto pz = dimuon.pt() * std::sinh(dimuon.eta());

      // muon variables
      auto py1 = dimuon.pt1() * ROOT::Math::sin(dimuon.phi1());
      auto px1 = dimuon.pt1() * ROOT::Math::cos(dimuon.phi1());
      auto py2 = dimuon.pt2() * ROOT::Math::sin(dimuon.phi2());
      auto px2 = dimuon.pt2() * ROOT::Math::cos(dimuon.phi2());

      // muon cuts

      // secondary vertexing test
      auto deltaZ = (dimuon.tauz() * pz) / dimuon.mass();

      if ((dimuon.eta1() > -4 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -4 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
        if (dimuon.pt1() >= ptCut && dimuon.pt2() >= ptCut) {
          if (dimuon.sign() == 0) {
            if (dimuon.chi2MatchMCHMFT1() <= 0. && dimuon.chi2MatchMCHMFT2() <= 0.) {
              if (rap > 2.5 && rap < 4) {
                // MCH-MID chi2 studies
                registry.get<TH1>(HIST("JPsiStandard"))->Fill(dimuon.mass());

                if (usePhi) {
                  if (dimuon.phi1() > PI / 2 || dimuon.phi1() < -PI / 2) {
                    hchi2MCHMIDLeft->Fill(dimuon.chi2MatchMCHMID1());
                  } else if (dimuon.phi1() < PI / 2 && dimuon.phi1() > -PI / 2) {
                    hchi2MCHMIDRight->Fill(dimuon.chi2MatchMCHMID1());
                  }

                  if ((dimuon.phi1() > PI / 2 || dimuon.phi1() < -PI / 2) && (dimuon.phi2() > PI / 2 || dimuon.phi2() < -PI / 2)) {
                    registry.get<TH1>(HIST("JPsiXNeg"))->Fill(dimuon.mass());
                  }

                  if ((dimuon.phi1() < PI / 2 && dimuon.phi1() > -PI / 2) && (dimuon.phi2() < PI / 2 && dimuon.phi2() > -PI / 2)) {
                    registry.get<TH1>(HIST("JPsiXPos"))->Fill(dimuon.mass());
                  }

                  if (dimuon.phi1() > PI / 2 || dimuon.phi1() < -PI / 2) {
                    if (dimuon.phi2() < PI / 2 && dimuon.phi2() > -PI / 2) {
                      registry.get<TH1>(HIST("JPsiXNeutral"))->Fill(dimuon.mass());
                    }
                  }
                  if (dimuon.phi1() < PI / 2 && dimuon.phi1() > -PI / 2) {
                    if (dimuon.phi2() > PI / 2 || dimuon.phi2() < -PI / 2) {
                      registry.get<TH1>(HIST("JPsiXNeutral"))->Fill(dimuon.mass());
                    }
                  }
                } // end if phi method

                if (useMomentum) {
                  if (px1 < 0. && px2 < 0.) {
                    hchi2MCHMIDLeft->Fill(dimuon.chi2MatchMCHMID1());
                  } else if (px1 > 0. && px2 > 0.) {
                    hchi2MCHMIDRight->Fill(dimuon.chi2MatchMCHMID1());
                  }
                  if (py1 < 0. && py2 < 0.) {
                    hchi2MCHMIDDown->Fill(dimuon.chi2MatchMCHMID1());
                  } else if (py1 > 0. && py2 > 0.) {
                    hchi2MCHMIDUp->Fill(dimuon.chi2MatchMCHMID1());
                  }
                  if (px1 < 0. && px2 < 0.) {
                    registry.get<TH1>(HIST("JPsiXNeg"))->Fill(dimuon.mass());
                  }
                  if (px1 > 0. && px2 > 0.) {
                    registry.get<TH1>(HIST("JPsiXPos"))->Fill(dimuon.mass());
                  }
                  if ((px1 < 0. && px2 > 0.) || (px1 > 0. && px2 < 0.)) {
                    registry.get<TH1>(HIST("JPsiXNeutral"))->Fill(dimuon.mass());
                  }
                  if (py1 < 0. && py2 < 0.) {
                    registry.get<TH1>(HIST("JPsiYNeg"))->Fill(dimuon.mass());
                  }
                  if (py1 > 0. && py2 > 0.) {
                    registry.get<TH1>(HIST("JPsiYPos"))->Fill(dimuon.mass());
                  }
                  if ((py1 < 0. && py2 > 0.) || (py1 > 0. && py2 < 0.)) {
                    registry.get<TH1>(HIST("JPsiYNeutral"))->Fill(dimuon.mass());
                  }
                } // end of momentum method

              } // end rapidity cut
            }
          } // end MCH-MID track selection
        } //end pt cut   
        // end OS selection
      }     // end eta cut (MUON acceptance)

      // MCH-MFT chi2 studies
      if ((dimuon.eta1() > -3.6 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -3.6 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
        if (dimuon.sign() == 0) {
          hdeltaZ->Fill(deltaZ);
        }

        if (dimuon.chi2MatchMCHMFT1() > 0 && dimuon.chi2MatchMCHMFT2() > 0) {
          hRapidityMCHMFT->Fill(rap);
          if (dimuon.sign() == 0) {
            if (rap > 2.5 && rap < 3.6) {
              // TODO check : MCH-MFT matching chi2 distribution for left and right MCH part
              if (dimuon.chi2MatchMCHMFT1() < 45 && dimuon.chi2MatchMCHMFT2() < 45) {
                hTauz->Fill(dimuon.tauz());
              }
              registry.get<TH3>(HIST("massTH3"))->Fill(dimuon.mass(), dimuon.chi2MatchMCHMFT1(), dimuon.chi2MatchMCHMFT2());
            } // end rapidity cut
          }   // end sign selection
        }     // end MCH-MFT chi2 selection
      }       // end eta cut

      // check on distribution without cuts
      hChi2MCHMFT->Fill(dimuon.chi2MatchMCHMFT1());
      hRapidityAll->Fill(rap);
    } // end loop over dimuons
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mumuReader>(cfgc)};
}
