// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFXToJpsiPiPiCandidateSelector.cxx
/// \brief X(3872) selection task.
/// \note Adapted from HFJpsiToEECandidateSelector.cxx
/// \author Rik Spijkers <r.spijkers@students.uu.nl>, Utrecht University

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "AnalysisCore/HFSelectorCuts.h"
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"
using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_x;
using namespace o2::analysis;

static const int npTBins = 9;
static const int nCutVars = 5;
//TODO: move this to namespace o2::analysis::hf_cuts_jpsi_toee, like done with the Jpsi selector
//    mass  dcaxy dcaz pt_jpsi pt_pi
constexpr double cuts[npTBins][nCutVars] =
  {{0.5, 0.2, 0.4, 1, 0},  /* pt<0.5   */
   {0.5, 0.2, 0.4, 1, 0},  /* 0.5<pt<1 */
   {0.5, 0.2, 0.4, 1, 0},  /* 1<pt<2   */
   {0.5, 0.2, 0.4, 1, 0},  /* 2<pt<3   */
   {0.5, 0.2, 0.4, 1, 0},  /* 3<pt<4   */
   {0.5, 0.2, 0.4, 1, 0},  /* 4<pt<5   */
   {0.5, 0.2, 0.4, 1, 0},  /* 5<pt<7   */
   {0.5, 0.2, 0.4, 1, 0},  /* 7<pt<10  */
   {0.5, 0.2, 0.4, 1, 0}}; /* 10<pt<15 */

/// Struct for applying Jpsi selection cuts
struct HFXToJpsiPiPiCandidateSelector {

  Produces<aod::HFSelXToJpsiPiPiCandidate> hfSelXToJpsiPiPiCandidate;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 50., "Upper bound of candidate pT"};

  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 10., "Upper bound of track pT for TPC PID"};

  Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};

  /// Gets corresponding pT bin from cut file array
  /// \param candpT is the pT of the candidate
  /// \return corresponding bin number of array
  template <typename T>
  int getpTBin(T candpT)
  {
    double pTBins[npTBins + 1] = {0, 0.5, 1., 2., 3., 4., 5., 7., 10., 15.};
    if (candpT < pTBins[0] || candpT >= pTBins[npTBins]) {
      return -1;
    }
    for (int i = 0; i < npTBins; i++) {
      if (candpT < pTBins[i + 1]) {
        return i;
      }
    }
    return -1;
  }

  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& track)
  {
    return true;
  }

  /// Conjugate independent toplogical cuts
  /// \param hfCandX is candidate
  /// \param trackNeutral is the track with the pi+ hypothesis
  /// \param trackPos is the track with the pi+ hypothesis
  /// \param trackNeg is the track with the pi- hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandX, const T2& hfCandJpsi, const T3& trackPos, const T3& trackNeg)
  {
    auto candpT = hfCandX.pt();
    int pTBin = getpTBin(candpT);
    if (pTBin == -1) {
      return false;
    }

    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false; //check that the candidate pT is within the analysis range
    }

    // TODO: replace hardcoded mass with "RecoDecay::getMassPDG(9920443)"
    if (TMath::Abs(InvMassXToJpsiPiPi(hfCandX) - 3872.) > cuts[pTBin][0]) {
      return false; //check that mass difference is within bounds
    }

    if ((hfCandJpsi.pt() < cuts[pTBin][3]) || (trackNeg.pt() < cuts[pTBin][4]) || (trackPos.pt() < cuts[pTBin][4])) {
      return false; //cut on daughter pT
    }
    // FIXME: dcaPrim0, dcaPrim1 are not defined???
    if (TMath::Abs(trackNeg.dcaPrim0()) > cuts[pTBin][1] || TMath::Abs(trackNeg.dcaPrim0()) > cuts[pTBin][1] || TMath::Abs(trackPos.dcaPrim0()) > cuts[pTBin][1]) {
      return false; //cut on daughter dca - need to add secondary vertex constraint here
    }
    if (TMath::Abs(trackNeg.dcaPrim1()) > cuts[pTBin][2] || TMath::Abs(trackPos.dcaPrim1()) > cuts[pTBin][2]) {
      return false; //cut on daughter dca - need to add secondary vertex constraint here
    }

    return true;
  }

  /// Check if track is ok for TPC PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TPC PID
  template <typename T>
  bool validTPCPID(const T& track)
  {
    if (TMath::Abs(track.pt()) < d_pidTPCMinpT || TMath::Abs(track.pt()) >= d_pidTPCMaxpT) {
      return false;
    }
    //if (track.TPCNClsFindable() < d_TPCNClsFindablePIDCut) return false;
    return true;
  }

  /// Check if track is compatible with given TPC Nsigma cut for the pion hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nsigma threshold to test against
  /// \return true if track satisfies TPC PID hypothesis for given Nsigma cut
  template <typename T>
  bool selectionPIDTPC(const T& track, int nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tpcNSigmaEl() < nSigmaCut;
  }
  /// PID selection on daughter track
  /// \param track is the daughter track
  /// \return 1 if successful PID match, 0 if successful PID rejection, -1 if no PID info
  template <typename T>
  int selectionPID(const T& track)
  {

    if (validTPCPID(track)) {
      if (!selectionPIDTPC(track, d_nSigmaTPC)) {

        return 0; //rejected by PID
      } else {
        return 1; //positive PID
      }
    } else {
      return -1; //no PID info
    }
  }
  void process(aod::hfCandX const& hfCandXs, aod::BigTracksPID const& tracks)
  {
    for (auto& hfCandX : hfCandXs) { //looping over X candidates
      // note the difference between Jpsi (index0) and pions (index1,2)
      auto candJpsi = hfCandX.index0(); // not a track
      auto trackPos = hfCandX.index1_as<aod::BigTracksPID>(); //positive daughter
      auto trackNeg = hfCandX.index2_as<aod::BigTracksPID>(); //negative daughter

      // check if flagged as X --> Jpsi Pi Pi
      if (!(hfCandX.hfflag() & 1 << XToJpsiPiPi)) {
        hfSelXToJpsiPiPiCandidate(0);
        continue;
      }

      // daughter track validity selection
      if (!daughterSelection(trackPos) || !daughterSelection(trackNeg)) {
        hfSelXToJpsiPiPiCandidate(0);
        continue;
      }

      //implement filter bit 4 cut - should be done before this task at the track selection level
      //need to add special cuts (additional cuts on decay length and d0 norm)

      if (!selectionTopol(hfCandX, trackPos, trackNeg)) {
        hfSelXToJpsiPiPiCandidate(0);
        continue;
      }

      if (selectionPID(trackPos) == 0 || selectionPID(trackNeg) == 0) {
        hfSelXToJpsiPiPiCandidate(0);
        continue;
      }

      hfSelXToJpsiPiPiCandidate(1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFXToJpsiPiPiCandidateSelector>(cfgc, TaskName{"hf-x-tojpsipipi-candidate-selector"})};
}
