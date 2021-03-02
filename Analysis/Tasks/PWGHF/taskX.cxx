// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskX.cxx
/// \brief X(3872) analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Rik Spijkers <r.spijkers@students.uu.nl>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"
#include "AnalysisCore/trackUtilities.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "ReconstructionDataFormats/DCA.h"
#include "AnalysisDataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/V0.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

namespace o2::aod
{
namespace extra
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
}
DECLARE_SOA_TABLE(Colls, "AOD", "COLLSID", o2::aod::extra::CollisionId);
} // namespace o2::aod
struct AddCollisionId {
  Produces<o2::aod::Colls> colls;
  void process(aod::HfCandProng2 const& candidates, aod::Tracks const&)
  {
    for (auto& candidate : candidates) {
      colls(candidate.index0_as<aod::Tracks>().collisionId());
    }
  }
};

/// X analysis task
/// FIXME: Still need to remove track duplication!!!
struct TaskX {
  HistogramRegistry registry{
    "registry",
    {{"hmassJpsi", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}}},
     {"hmassX", "3-prong candidates;inv. mass (#J/psi pi+ pi-) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 2., 7.}}}},
     {"hptcand", "X candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}}}};

  Configurable<int> d_selectionFlagJpsi{"d_selectionFlagJpsi", 1, "Selection Flag for Jpsi"};
  Configurable<double> cutEtaCandMax{"cutEtaCandMax", -1., "max. cand. pseudorapidity"};

  // vertexing parameters
  Configurable<double> d_bz{"d_bz", 5., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_jpsi::isSelJpsiToEE >= d_selectionFlagJpsi);

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massJpsi = RecoDecay::getMassPDG(443);

  /// aod::BigTracks is not soa::Filtered, should be added when filters are added
  void process(aod::Collision const&, aod::BigTracks const& tracks, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelJpsiToEECandidate, aod::Colls>> const& candidates)
  {
    //Initialise fitter for X vertex
    o2::vertexing::DCAFitterN<3> Xfitter;
    Xfitter.setBz(d_bz);
    Xfitter.setPropagateToPCA(b_propdca);
    Xfitter.setMaxR(d_maxr);
    Xfitter.setMinParamChange(d_minparamchange);
    Xfitter.setMinRelChi2Change(d_minrelchi2change);
    Xfitter.setUseAbsDCA(d_UseAbsDCA);

    //Initial fitter to redo J/psi-vertex to get extrapolated daughter tracks
    o2::vertexing::DCAFitterN<2> Jpsifitter;
    Jpsifitter.setBz(d_bz);
    Jpsifitter.setPropagateToPCA(b_propdca);
    Jpsifitter.setMaxR(d_maxr);
    Jpsifitter.setMinParamChange(d_minparamchange);
    Jpsifitter.setMinRelChi2Change(d_minrelchi2change);
    Jpsifitter.setUseAbsDCA(d_UseAbsDCA);

    // loop over J/psi candidates (picked up from candidate selector)
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << JpsiToEE)) {
        continue;
      }
      if (cutEtaCandMax >= 0. && std::abs(candidate.eta()) > cutEtaCandMax) {
        continue;
      }
      registry.fill(HIST("hmassJpsi"), InvMassJpsiToEE(candidate));

      const std::array<float, 3> vertexJpsi = {candidate.xSecondaryVertex(), candidate.ySecondaryVertex(), candidate.zSecondaryVertex()};
      const std::array<float, 3> momentumJpsi = {candidate.px(), candidate.py(), candidate.pz()};

      // auto prong0 = candidate.index0(); //index0_as<aod::BigTracks>();
      // auto prong1 = candidate.index1(); //index1_as<aod::BigTracks>();
      // auto prong0TrackParCov = getTrackParCov(prong0);
      // auto prong1TrackParCov = getTrackParCov(prong1);

      int index0jpsi = candidate.index0Id();
      int index1jpsi = candidate.index1Id();
      auto JpsiProng0 = candidate.index0();
      auto JpsiProng1 = candidate.index1();
      auto prong0TrackParCov = getTrackParCov(JpsiProng0);
      auto prong1TrackParCov = getTrackParCov(JpsiProng1);

      // reconstruct Jpsi secondary vertex
      if (Jpsifitter.process(prong0TrackParCov, prong1TrackParCov) == 0) {
        continue;
      }

      //Propogate prong tracks to extrapolated track position.
      prong0TrackParCov.propagateTo(Jpsifitter.getTrack(0).getX(), d_bz);
      prong1TrackParCov.propagateTo(Jpsifitter.getTrack(1).getX(), d_bz);

      // build a Jpsi neutral track [Q]
      auto trackJpsi = o2::dataformats::V0(vertexJpsi, momentumJpsi, prong0TrackParCov, prong1TrackParCov, {0, 0}, {0, 0});

      for (auto PionTrackPos = tracks.begin(); PionTrackPos != tracks.end(); ++PionTrackPos) {
        if (PionTrackPos.signed1Pt() < 0) {
          continue;
        }
        if (PionTrackPos.globalIndex() == index0jpsi) {
          printf("pos track id check: %ld, %u\n", PionTrackPos.globalIndex(), index0jpsi);
          printf("pos track pt check: %f, %f\n", JpsiProng0.pt(), PionTrackPos.pt());
          continue;
        }
        for (auto PionTrackNeg = tracks.begin(); PionTrackNeg != tracks.end(); ++PionTrackNeg) {
          if (PionTrackNeg.signed1Pt() > 0) {
            continue;
          }
          if (PionTrackNeg.globalIndex() == index1jpsi) {
            printf("neg track id check: %ld, %u\n", PionTrackNeg.globalIndex(), index1jpsi);
            printf("neg track pt check: %f, %f\n", JpsiProng1.pt(), PionTrackNeg.pt());
            continue;
          }
          registry.fill(HIST("hptcand"), candidate.pt() + PionTrackPos.pt() + PionTrackNeg.pt());

          auto piplusTrackParCov = getTrackParCov(PionTrackPos);
          auto piminusTrackParCov = getTrackParCov(PionTrackNeg);

          std::array<float, 3> pvecJpsi = {0., 0., 0.};
          std::array<float, 3> pvecpiplus = {0., 0., 0.};
          std::array<float, 3> pvecpiminus = {0., 0., 0.};
          std::array<float, 3> pvecXCand = {0., 0., 0.};

          //find the DCA between the Jpsi and the bachelor track, for X
          int nCand = Xfitter.process(trackJpsi, piplusTrackParCov, piminusTrackParCov); //Plot nCand
          if (nCand == 0)
            continue;

          Xfitter.propagateTracksToVertex();             // propagate the pions and Jpsi to the X vertex
          Xfitter.getTrack(0).getPxPyPzGlo(pvecJpsi);    // momentum of Jpsi at the X vertex; Need to fill new Jpsi pT
          Xfitter.getTrack(1).getPxPyPzGlo(pvecpiplus);  // momentum of pi+ at the X vertex; Need to fill new p+ pT
          Xfitter.getTrack(2).getPxPyPzGlo(pvecpiminus); // momentum of pi- at the X vertex; Need to fill new p- pT

          pvecXCand = array{
            pvecJpsi[0] + pvecpiplus[0] + pvecpiminus[0],
            pvecJpsi[1] + pvecpiplus[1] + pvecpiminus[1],
            pvecJpsi[2] + pvecpiplus[2] + pvecpiminus[2],
          };

          // invariant mass
          auto arrMom = array{pvecJpsi, pvecpiplus, pvecpiminus};
          auto arrMass = array{massJpsi, massPi, massPi};
          double massX = RecoDecay::M(arrMom, arrMass); //Plot X mass
          registry.fill(HIST("hmassX"), massX);
        } // pi- loop
      }   // pi+ loop
    }     // Jpsi loop
  }       // process
};        // struct

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<AddCollisionId>("hf-task-add-collisionId"),
    adaptAnalysisTask<TaskX>("hf-task-x")};
  return workflow;
}
