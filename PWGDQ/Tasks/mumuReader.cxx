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
// #include "meta/inc/TStreamerElement.h"
#include <TH1F.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;

struct mumuReader {

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> chi2MCHMIDLeft{TH1F("chi2MCHMIDLeft", "chi2MCHMIDLeft", 100, 0., 16.)};
  OutputObj<TH1F> chi2MCHMIDRight{TH1F("chi2MCHMIDRight", "chi2MCHMIDRight", 100, 0., 16.)};
  OutputObj<TH1F> hRapidityAll{TH1F("Rapidity", "Rapidity", 267, 2., 4.)};
  OutputObj<TH1F> hRapidityMCHMFT{TH1F("RapidityMCHMFT", "RapidityMCHMFT", 267, 2., 4.)};

  // histogram registry test
  HistogramRegistry registry{
    "registry",
    {{"mass", "M", {HistType::kTH1F, {{750, 0, 15}}}}
     }};
  // end of histogram registry test

  void init(o2::framework::InitContext&)
  {
    AxisSpec massAxis = {750, 0, 15, "M"};
    AxisSpec chi2Axis1 = {200,0,200,"#chi^2_{MCH-MFT} for mu 1"};
    AxisSpec chi2Axis2 = {200,0,200,"#chi^2_{MCH-MFT} for mu 2"};
    HistogramConfigSpec histospec({HistType::kTH3F, {massAxis,chi2Axis1,chi2Axis2}});

    registry.add("massTH3","Mass TH3 Histogram", histospec);
  };

  void process(aod::DimuonsAll const& dimuons)
  {

    LOG(info) << "NEW FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    for (auto& dimuon : dimuons) {
      auto rap = -ROOT::Math::log((ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt() * ROOT::Math::cosh(dimuon.eta()) * ROOT::Math::cosh(dimuon.eta())) + dimuon.pt() * ROOT::Math::sinh(dimuon.eta())) / (ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt())));
      hRapidityAll->Fill(rap);

      if (dimuon.phi1() > 3.14 / 2 || dimuon.phi1() < -3.14 / 2) {
        chi2MCHMIDLeft->Fill(dimuon.chi2MatchMCHMID1());
      } else if (dimuon.phi1() < 3.14 / 2 && dimuon.phi1() > -3.14 / 2) {
        chi2MCHMIDRight->Fill(dimuon.chi2MatchMCHMID1());
      }                                                                                                       // end MCH MID chi2 check
      if ((dimuon.eta1() > -3.6 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -3.6 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
        if (dimuon.chi2MatchMCHMFT1() > 0 && dimuon.chi2MatchMCHMFT2() > 0) {
          hRapidityMCHMFT->Fill(rap);
          if (dimuon.sign() == 0) {
            registry.get<TH3>(HIST("massTH3"))->Fill(dimuon.mass(),dimuon.chi2MatchMCHMFT1(),dimuon.chi2MatchMCHMFT2());
          } // end sign selection
        }   // end MCH-MFT chi2 selection
      }     // end eta cut
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mumuReader>(cfgc)};
}
