#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelectedPrefilter, isMuonSelectedPrefilter, int);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);//, dqanalysisflags::IsMuonSelectedPrefilter);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsPrefilterVetoed);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;

constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;

constexpr static int pairTypeMuMu = VarManager::kDecayToMuMu;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  // Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  // Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};

  HistogramManager* fHistMan = nullptr;
  // MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;

  void init(o2::framework::InitContext&)
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    VarManager::SetDefaultVarNames();
    /*
    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;", fConfigAddEventHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                           // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }*/
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent const& event)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<TEventFillMap>(event);
    // TODO: make this condition at compile time
    /*if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
    }*/
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      /*
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }
      */
      eventSel(1);
    } else {
      eventSel(0);
    }
  }

  void processSkimmed(MyEvents::iterator const& event)
  {
    runEventSelection<gkEventFillMap>(event);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", false);
  // TODO: Add process functions subscribing to Framework Collision
};

struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fMuonCuts;

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    /*if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackMuon_BeforeCuts;";
      for (auto& cut : fMuonCuts) {
        histDirNames += Form("TrackMuon_%s;", cut.GetName());
      }

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddMuonHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                          // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }*/
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runMuonSelection(TEvent const& event, TMuons const& muons)
  {
    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
    VarManager::FillEvent<TEventFillMap>(event);

    muonSel.reserve(muons.size());
    uint32_t filterMap = 0;
    int iCut = 0;

    for (auto& muon : muons) {
      filterMap = 0;
      VarManager::FillTrack<TMuonFillMap>(muon);
      /*if (fConfigQA) { // TODO: make this compile time
        fHistMan->FillHistClass("TrackMuon_BeforeCuts", VarManager::fgValues);
      }*/

      iCut = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << iCut);
          /*if (fConfigQA) { // TODO: make this compile time
            fHistMan->FillHistClass(Form("TrackMuon_%s", (*cut).GetName()), VarManager::fgValues);
          }*/
        }
      }
      muonSel(static_cast<int>(filterMap));
    } // end loop over tracks
  }

  void processSkimmed(MyEvents::iterator const& event, MyMuonTracks const& muons)
  {
    runMuonSelection<gkEventFillMap, gkMuonFillMap>(event, muons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", false);
};

struct AnalysisSameEventPairing {

  Produces<aod::DimuonsAll> dimuonAllList;
  //Produces<aod::MuonList> MuonList;
  Produces<aod::DileptonFlow> dileptonFlowList;
  float mMagField = 0.0;
  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};
  Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<bool> fConfigFlatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected == 1;
  // NOTE: the barrel filter map contains decisions for both electrons and hadrons used in the correlation task
  Filter filterMuonTrackSelected = aod::dqanalysisflags::isMuonSelected > 0;
  //Filter prefilter = aod::dqanalysisflags::isPrefilterVetoed == 0;

  HistogramManager* fHistMan;

  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  // TODO: provide as Configurable the list and names of the cuts which should be used in pairing
  uint32_t fTwoTrackFilterMask = 0;
  uint32_t fTwoMuonFilterMask = 0;
  std::vector<std::vector<TString>> fMuonHistNames;
  std::vector<AnalysisCompositeCut> fPairCuts;

  void init(o2::framework::InitContext& context)
  {
    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Keep track of all the histogram class names to avoid composing strings in the event mixing pairing
    TString histNames = "";
    std::vector<TString> names;

    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    if (context.mOptions.get<bool>("processDecayToMuMuSkimmed") || context.mOptions.get<bool>("processDecayToMuMuVertexingSkimmed")) {
      TString cutNames = fConfigMuonCuts.value;
      if (!cutNames.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
        for (int icut = 0; icut < objArray->GetEntries(); ++icut) { // loop over track cuts
          fTwoMuonFilterMask |= (uint32_t(1) << icut);
          // no pair cuts
          names = {
            Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
            Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
          histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
          fMuonHistNames.push_back(names);

          TString cutNamesStr = fConfigPairCuts.value;
          if (!cutNamesStr.IsNull()) { // if pair cuts
            std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
            for (int iPairCut = 0; iPairCut < objArrayPair->GetEntries(); ++iPairCut) { // loop over pair cuts
              names = {
                Form("PairsMuonSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                Form("PairsMuonSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                Form("PairsMuonSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
              histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
              fMuonHistNames.push_back(names);
            } // end loop (pair cuts)
          }   // end if (pair cuts)
        }     // end loop (muon cuts)
      }       // end if (muon cuts)
    }

    // Usage example of ccdb
    // ccdb->setURL(url.value);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setCreatedNotAfter(nolaterthan.value);

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                      // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks1, typename TTracks2>
  void runSameEventPairing(TEvent const& event, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if (fCurrentRun != event.runNumber()) {
      if (fUseRemoteField.value) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, event.timestamp());
        if (grpmag != nullptr) {
          mMagField = grpmag->getNominalL3Field();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", event.timestamp());
        }
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(mMagField);
        } else {
          VarManager::SetupTwoProngFwdDCAFitter(mMagField, fPropToPCA.value, 200.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
        }
      } else {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupTwoProngFwdDCAFitter(fConfigMagField.value, fPropToPCA.value, 200.0f, 1.0e-3f, 0.9f, fUseAbsDCA.value);
        }
      }
      fCurrentRun = event.runNumber();
    }

    TString cutNames = fConfigMuonCuts.value;
    std::vector<std::vector<TString>> histNames = fMuonHistNames;

    std::unique_ptr<TObjArray> objArray(cutNames.Tokenize(","));
    int ncuts = objArray->GetEntries();

    uint32_t twoTrackFilter = 0;
    uint32_t dileptonFilterMap = 0;
    if (fConfigFlatTables.value) {
      dimuonAllList.reserve(1);
      //MuonList.reserve(1);
    }
    for (auto& [t1, t2] : combinations(tracks1, tracks2)) {
      if constexpr (TPairType == VarManager::kDecayToMuMu) {
        twoTrackFilter = uint32_t(t1.isMuonSelected()) & uint32_t(t2.isMuonSelected()) & fTwoMuonFilterMask;
      }
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }

      // TODO: FillPair functions need to provide a template argument to discriminate between cases when cov matrix is available or not
      VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
      if constexpr (TPairType == pairTypeMuMu) { // call this just for ee or mumu pairs
        VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2);
      }

      // TODO: provide the type of pair to the dilepton table (e.g. ee, mumu, emu...)
      dileptonFilterMap = twoTrackFilter;

      constexpr bool muonHasCov = ((TTrackFillMap & VarManager::ObjTypes::MuonCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedMuonCov) > 0);
      if constexpr ((TPairType == pairTypeMuMu) && muonHasCov) {
        if (fConfigFlatTables.value) {
          dimuonAllList(event.posX(), event.posY(), event.posZ(), -999., -999., -999., VarManager::fgValues[VarManager::kMass], false, VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], t1.pt(), t1.eta(), t1.phi(), t1.sign(), t2.pt(), t2.eta(), t2.phi(), t2.sign(), 0., 0., t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(), t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(), -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., t1.isAmbiguous(), t2.isAmbiguous());
          //MuonList(VarManager::fgValues[VarManager::kMass],t1.sign()+t2.sign(),t1.x(),t2.x(),t1.y(),t2.y(),t1.z(),t2.z(),t1.pt(), t1.eta(), t1.phi(), t1.sign(), t2.pt(), t2.eta(), t2.phi(), t2.sign(), t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(), t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT())
        }
      }

      int iCut = 0;
      for (int icut = 0; icut < ncuts; icut++) {
        if (twoTrackFilter & (uint32_t(1) << icut)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(histNames[iCut][0].Data(), VarManager::fgValues);
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(histNames[iCut][1].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(histNames[iCut][2].Data(), VarManager::fgValues);
            }
          }
          iCut++;
          for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++, iCut++) {
            AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
            if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
              continue;
            if (t1.sign() * t2.sign() < 0) {
              fHistMan->FillHistClass(histNames[iCut][0].Data(), VarManager::fgValues);
            } else {
              if (t1.sign() > 0) {
                fHistMan->FillHistClass(histNames[iCut][1].Data(), VarManager::fgValues);
              } else {
                fHistMan->FillHistClass(histNames[iCut][2].Data(), VarManager::fgValues);
              }
            }
          }      // end loop (pair cuts)
        } else { // end if (filter bits)
          iCut++;
        }
      } // end loop (cuts)
    }   // end loop over pairs
  }

  void processDecayToMuMuSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMap>(event, muons, muons);
  }
  void processDecayToMuMuVertexingSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyMuonTracksSelectedWithCov> const& muons)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
    runSameEventPairing<VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, muons, muons);
  }

  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuSkimmed, "Run muon-muon pairing, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDecayToMuMuVertexingSkimmed, "Run muon-muon pairing and vertexing, with skimmed muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", false);
};


struct SecondVertexing {
  SliceCache cache;
  Produces<aod::Prefilter> prefilter;
  Preslice<MyMuonTracks> perCollision = aod::reducedtrack::reducedeventId;

  // Configurables
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};

  Filter muonTracksSelectedPrefilter = aod::dqanalysisflags::isMuonSelectedPrefilter > 0;

  Partition<soa::Filtered<MyMuonTracksSelected>> muonTracksSelected = aod::dqanalysisflags::isMuonSelected > 0;

  std::map<int, bool> fPrefiltermap;
  AnalysisCompositeCut* fPairCut;

  void init(o2::framework::InitContext& context)
  {
    cout << "Hi !" << endl;
    fPairCut = new AnalysisCompositeCut(true);
    TString pairCutStr = fConfigPrefilterPairCut.value;
    if (!pairCutStr.IsNull()) {
      fPairCut = dqcuts::GetCompositeCut(pairCutStr.Data());
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <int TPairType, uint32_t TTrackFillMap, typename TTracks1, typename TTracks2>
  void runPrefilterPairing(TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      if (track1.sign() * track2.sign() > 0) {
        continue;
      }
      cout << "Hi !" << endl;

      // pairing
      VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);

      if (fPairCut->IsSelected(VarManager::fgValues)) {
        fPrefiltermap[track1.globalIndex()] = true;
        fPrefiltermap[track2.globalIndex()] = true;
      }
    }
  }

  void processMuonSkimmed(MyEventsSelected const& events, soa::Filtered<MyMuonTracksSelected> const& filteredTracks, MyMuonTracks const& tracks)
  {
    // cout << "Hi !" << endl;
    const int pairType = VarManager::kDecayToMuMu;
    fPrefiltermap.clear();

    for (auto& event : events) {
      if (event.isEventSelected()) {
        auto groupedPrefilterCandidates = filteredTracks.sliceBy(perCollision, event.globalIndex());
        auto groupedMuonCandidates = muonTracksSelected->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex(), cache);
        runPrefilterPairing<pairType, gkMuonFillMap>(groupedPrefilterCandidates, groupedMuonCandidates);
      }
    } // end loop events

    // Fill Prefilter bits for all tracks to have something joinable to MyBarrelTracksSelected
    for (auto& track : tracks) {
      prefilter(static_cast<int>(fPrefiltermap[track.globalIndex()]));
    }
  }

  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(SecondVertexing, processMuonSkimmed, "Run Prefilter selection on reduced tracks", false);
  PROCESS_SWITCH(SecondVertexing, processDummy, "Do nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc)};
    //adaptAnalysisTask<SecondVertexing>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    TString histName = configVar.value;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
        if (classStr.Contains("PIDCalibElectron")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PIDCalibPion")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PIDCalibProton")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-correlation");
    }
  } // end loop over histogram classes
}