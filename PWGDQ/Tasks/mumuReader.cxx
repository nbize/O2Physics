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
  OutputObj<TH1F> hChi2MCHMFT{TH1F("Chi2MCHMFT", "Chi2MCHMFT", 100, 0., 200.)};
  OutputObj<TH1F> hchi2MCHMFTLeft{TH1F("Chi2MCHMFTLeft", "Chi2MCHMFTLeft", 200, 0., 200.)};
  OutputObj<TH1F> hchi2MCHMFTRight{TH1F("Chi2MCHMFTRight", "Chi2MCHMFTRight", 200, 0., 200.)};

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
    AxisSpec ptAxis = {120, 0.0, 30.0, "p_{T}"};
    AxisSpec rapAxis = {200, 2.5, 4.0, "y"};
    AxisSpec ambiAxis = {50, -1, 4, "isAmbi"};
    AxisSpec chi2Axis1 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 1"};
    AxisSpec chi2Axis2 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 2"};
    AxisSpec trackChi2Axis = {200, 0.0, 200.0, "#chi^2"};
    HistogramConfigSpec histospec({HistType::kTH3F, {massAxis, chi2Axis1, chi2Axis2}});
    HistogramConfigSpec massSpec({HistType::kTH1F, {massAxis}});
    HistogramConfigSpec ptSpec({HistType::kTH1F, {ptAxis}});
    HistogramConfigSpec ambiSpec({HistType::kTH1F, {ambiAxis}});
    HistogramConfigSpec trackChi2Spec({HistType::kTH1F, {trackChi2Axis}});
    HistogramConfigSpec rapSpec({HistType::kTH1F, {rapAxis}});
    HistogramConfigSpec massChi2Spec({HistType::kTH2F, {massAxis, trackChi2Axis}});
    // TODO : create THnSparse to avoid to put the cuts directly in the task
    HistogramConfigSpec sparseSpec({HistType::kTHnSparseF, {massAxis, ambiAxis, ambiAxis, chi2Axis1, chi2Axis2}});

    // MUON standalone studies
    registry.add("JPsiXPos", "JPsi X pos", massSpec);
    registry.add("JPsiXNeg", "JPsi X neg", massSpec);
    registry.add("JPsiXNeutral", "JPsi X neutral", massSpec);
    registry.add("JPsiYPos", "JPsi Y pos", massSpec);
    registry.add("JPsiYNeg", "JPsi Y neg", massSpec);
    registry.add("JPsiYNeutral", "JPsi Y neutral", massSpec);
    registry.add("JPsiStandard", "JPsi standard", massSpec);
    registry.add("JPsiStandardNoAmbig", "JPsi standard without ambiguous tracks", massSpec);
    registry.add("JPsiStandardAmbigOnly", "JPsi standard with ambiguous tracks only", massSpec);
    registry.add("TH2MassChi2", "Mass X chi2", massChi2Spec);

    // pt bins
    registry.add("Pt", "pt distribution", ptSpec);

    registry.add("JPsiPt01", "JPsi pt 0 - 1", massSpec);
    registry.add("JPsiPt12", "JPsi pt 1 - 2", massSpec);
    registry.add("JPsiPt23", "JPsi pt 2 - 3", massSpec);

    registry.add("JPsiPt34", "JPsi pt 3 - 4", massSpec);
    registry.add("JPsiPt45", "JPsi pt 4 - 5", massSpec);
    registry.add("JPsiPt58", "JPsi pt 5 - 8", massSpec);
    registry.add("JPsiPt812", "JPsi pt 8 - 12", massSpec);
    registry.add("JPsiPt820", "JPsi pt 8 - 20", massSpec);

    registry.add("JPsiPt35", "JPsi pt 3 - 5", massSpec);
    registry.add("JPsiPt57", "JPsi pt 5 - 7", massSpec);
    registry.add("JPsiPt710", "JPsi pt 7 - 10", massSpec);
    registry.add("JPsiPt1020", "JPsi pt 10 - 20", massSpec);

    // rapidity bins
    registry.add("JPsiRap253", "JPsi y 2.5 - 3", massSpec);
    registry.add("JPsiRap335", "JPsi y 3 - 3.5", massSpec);
    registry.add("JPsiRap354", "JPsi y 3.5 - 4", massSpec);

    // ambiguous tracks
    registry.add("ambiguous1", "ambiguous mu 1", ambiSpec);
    registry.add("ambiguous2", "ambiguous mu 2", ambiSpec);

    // chi2 track
    registry.add("trackChi2_1", "chi2 of track distribution", trackChi2Spec);
    registry.add("trackChi2_2", "chi2 of track distribution", trackChi2Spec);
    registry.add("trackChi2_1_noAmbig", "chi2 of track distribution (no ambig)", trackChi2Spec);
    // MCH-MFT chi2 studies
    registry.add("massTH3", "Mass TH3 Histogram", histospec);
    registry.add("ambigGlobalJPsi", "THnSparse Jpsi", sparseSpec);
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

                //______________________________________________________________________________
                // track chi2 studies
                registry.get<TH1>(HIST("trackChi2_1"))->Fill(dimuon.chi21());
                registry.get<TH1>(HIST("trackChi2_2"))->Fill(dimuon.chi22());
                registry.get<TH2>(HIST("TH2MassChi2"))->Fill(dimuon.mass(), dimuon.chi21());

                //______________________________________________________________________________
                // ambiguous tracks studies
                registry.get<TH1>(HIST("ambiguous1"))->Fill(dimuon.isAmbig1());
                registry.get<TH1>(HIST("ambiguous2"))->Fill(dimuon.isAmbig2());
                if (!(dimuon.isAmbig1()) || !(dimuon.isAmbig2())) {
                  registry.get<TH1>(HIST("JPsiStandardNoAmbig"))->Fill(dimuon.mass());
                  registry.get<TH1>(HIST("trackChi2_1_noAmbig"))->Fill(dimuon.chi21());
                  //______________________________________________________________________________
                  // pt bins
                  if (dimuon.pt() > 0 && dimuon.pt() < 1) {
                    registry.get<TH1>(HIST("JPsiPt01"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 1 && dimuon.pt() < 2) {
                    registry.get<TH1>(HIST("JPsiPt12"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 2 && dimuon.pt() < 3) {
                    registry.get<TH1>(HIST("JPsiPt23"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 3 && dimuon.pt() < 4) {
                    registry.get<TH1>(HIST("JPsiPt34"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 4 && dimuon.pt() < 5) {
                    registry.get<TH1>(HIST("JPsiPt45"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 5 && dimuon.pt() < 8) {
                    registry.get<TH1>(HIST("JPsiPt58"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 8 && dimuon.pt() < 12) {
                    registry.get<TH1>(HIST("JPsiPt812"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 8 && dimuon.pt() < 20) {
                    registry.get<TH1>(HIST("JPsiPt820"))->Fill(dimuon.mass());
                  }
                  // different pt bins
                  if (dimuon.pt() > 3 && dimuon.pt() < 5) {
                    registry.get<TH1>(HIST("JPsiPt35"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 5 && dimuon.pt() < 7) {
                    registry.get<TH1>(HIST("JPsiPt57"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 7 && dimuon.pt() < 10) {
                    registry.get<TH1>(HIST("JPsiPt710"))->Fill(dimuon.mass());
                  }
                  if (dimuon.pt() > 10 && dimuon.pt() < 20) {
                    registry.get<TH1>(HIST("JPsiPt1020"))->Fill(dimuon.mass());
                  }

                  //______________________________________________________________________________
                  // rapidity bins
                  if (rap > 2.5 && rap < 3.) {
                    registry.get<TH1>(HIST("JPsiRap253"))->Fill(dimuon.mass());
                  }
                  if (rap > 3. && rap < 3.5) {
                    registry.get<TH1>(HIST("JPsiRap335"))->Fill(dimuon.mass());
                  }
                  if (rap > 3.5 && rap < 4.) {
                    registry.get<TH1>(HIST("JPsiRap354"))->Fill(dimuon.mass());
                  }
                }
                if (dimuon.isAmbig1() || dimuon.isAmbig2()) {
                  registry.get<TH1>(HIST("JPsiStandardAmbigOnly"))->Fill(dimuon.mass());
                }

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
        }   // end pt cut
        // end OS selection
      } // end eta cut (MUON acceptance)

      // MCH-MFT chi2 studies
      if ((dimuon.eta1() > -3.6 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -3.6 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
        if (dimuon.sign() == 0) {
          hdeltaZ->Fill(deltaZ);
        }

        if (dimuon.chi2MatchMCHMFT1() > 0 && dimuon.chi2MatchMCHMFT2() > 0) {
          if (!(dimuon.isAmbig1()) && !(dimuon.isAmbig2())) { //
            hRapidityMCHMFT->Fill(rap);
            if (dimuon.sign() == 0) {
              if (rap > 2.5 && rap < 3.6) {
                // TODO check : MCH-MFT matching chi2 distribution for left and right MCH part
                if (dimuon.chi2MatchMCHMFT1() < 45 && dimuon.chi2MatchMCHMFT2() < 45) {
                  hTauz->Fill(dimuon.tauz());
                }
                if (px1 < 0. && px2 < 0.) {
                  hchi2MCHMFTLeft->Fill(dimuon.chi2MatchMCHMFT1());
                } else if (px1 > 0. && px2 > 0.) {
                  hchi2MCHMFTRight->Fill(dimuon.chi2MatchMCHMFT1());
                }
                registry.get<TH3>(HIST("massTH3"))->Fill(dimuon.mass(), dimuon.chi2MatchMCHMFT1(), dimuon.chi2MatchMCHMFT2());
                registry.get<THnSparse>(HIST("ambigGlobalJPsi"))->Fill(dimuon.mass(), dimuon.isAmbig1(), dimuon.isAmbig2(), dimuon.chi2MatchMCHMFT1(), dimuon.chi2MatchMCHMFT2());
              } // end rapidity cut
            }   // end sign selection
          }     // ambig cut
        }       // end MCH-MFT chi2 selection
      }         // end eta cut

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