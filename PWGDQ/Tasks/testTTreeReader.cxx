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

struct ReadTTree {
  // void init() {};
  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> chi2MCHMIDLeft{TH1F("chi2MCHMIDLeft", "chi2MCHMIDLeft", 100, 0., 16.)};
  OutputObj<TH1F> chi2MCHMIDRight{TH1F("chi2MCHMIDRight", "chi2MCHMIDRight", 100, 0., 16.)};
  void process(aod::DimuonsAll const& dimuons)
  {
    for (auto& dimuon : dimuons) {
      //LOG(info) << "This dimuon has a of pt = " << dimuon.pt() << "and an mass of m = " << dimuon.mass();
      if (dimuon.phi1()> 3.14/2 || dimuon.phi2() < -3.14/2){
        chi2MCHMIDLeft->Fill(dimuon.chi2MatchMCHMID2());
      }
      else if (dimuon.phi1() < 3.14/2 && dimuon.phi2() > -3.14/2){
        chi2MCHMIDRight->Fill(dimuon.chi2MatchMCHMID1());
      }
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ReadTTree>(cfgc),
  };
}
