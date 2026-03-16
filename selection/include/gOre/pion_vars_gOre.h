/**
 * @file pion_vars_gOre.h
 * @brief Truth variables for studying misidentified NC pi0 events (true_mc_category == 2).
 *
 * Provides three groups of variables:
 *
 *  1. MCTruth scope  — accesses the pi0 GENIE primary from obj.prim (SRTrueParticle).
 *     Full G4 info available: generation/start/end positions, wallout/wallin, process
 *     codes (as g4_process_ integers), generator, contained, vis_e, nhit, n_daughters.
 *
 *  2. Selectors (registered for both true/reco; only meaningful in true mode):
 *       leading_pi0_photon    — highest-KE SPINE truth photon with ancestor_pdg_code==111
 *       subleading_pi0_photon — second-highest-KE SPINE truth photon with same criterion
 *     Use in TOML as:  type = "true",  selector = "leading_pi0_photon"
 *     Branch names produced:  true_leading_pi0_photon_<varname>
 *                             true_subleading_pi0_photon_<varname>
 *
 *  3. TrueParticle scope — per-photon G4-truth variables from SRParticleTruthDLP.
 *     These are paired with the selectors above to produce one value per interaction.
 *     Available fields (compared to scan_pi0.C):
 *       pi0_photon_ke              ← pvars::ke(p)              [MeV]   (~startE - mass)
 *       pi0_photon_energy_init     ← p.energy_init             [MeV]   (total initial energy)
 *       pi0_photon_energy_deposit  ← p.energy_deposit          [MeV]   (like visE)
 *       pi0_photon_gen_x/y/z       ← p.position[0/1/2]         [cm]    (creation point = tp_gen)
 *       pi0_photon_start_x/y/z     ← p.start_point[0/1/2]      [cm]    (TPC entry  = tp_start)
 *       pi0_photon_end_x/y/z       ← p.end_point[0/1/2]        [cm]    (end point  = tp_end)
 *       pi0_photon_dir_x/y/z       ← p.start_dir[0/1/2]                (unit direction)
 *       pi0_photon_length          ← p.length                  [cm]
 *       pi0_photon_is_contained    ← p.is_contained                     (1/0)
 *       pi0_photon_is_primary      ← pvars::primary_classification(p)   (SPINE primary flag)
 *       pi0_photon_ancestor_pdg    ← p.ancestor_pdg_code                (111 for pi0 daughters)
 *       pi0_photon_track_id        ← p.track_id                         (G4 track ID)
 *       pi0_photon_parent_id       ← p.parent_id                        (parent G4 track ID)
 *       pi0_photon_px/py/pz        ← p.momentum[0/1/2]         [MeV/c]
 *     NOT available (SRTrueParticle only, not in SRParticleTruthDLP):
 *       wallout, wallin, start_process/end_process (int enum), generator (int enum)
 *       These are fully available for the pi0 ITSELF via the MCTruth scope above.
 *
 * Wall_t encoding used by pion_wallout / pion_wallin (MCTruth scope):
 *   0 = kWallNone, 1 = kWallTop, 2 = kWallBottom, 3 = kWallLeft,
 *   4 = kWallRight, 5 = kWallFront (upstream), 6 = kWallBack (downstream)
 *
 * g4_process_ encoding (pion_start_process, pion_end_process):
 *   0=kG4primary  3=kG4Decay  36=kG4conv  37=kG4phot  51=kG4compt  1024=kG4UNKNOWN
 *
 * generator_ encoding (pion_generator):
 *   0=kUnknownGenerator  1=kGENIE  2=kMeVPrtl
 **/

#ifndef PION_VARS_GORE_H
#define PION_VARS_GORE_H

#include "include/gOre/vars_gOre.h"

// ============================================================
//  SELECTORS — pi0 decay photons (true interactions only)
// ============================================================
namespace selectors::gOre
{
  /**
   * @brief Index of the leading (highest-KE) SPINE truth photon whose
   * ancestor is a pi0 (ancestor_pdg_code == 111).
   * Returns kNoMatch if none found. For reco interactions, always kNoMatch.
   **/
  template <class T>
  size_t leading_pi0_photon(const T& obj)
  {
    size_t best_idx = kNoMatch;
    double best_ke  = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < obj.particles.size(); ++i)
    {
      const auto& p = obj.particles[i];
      bool is_pi0_photon = false;
      if constexpr (std::is_same_v<T, TType>)
        is_pi0_photon = (pvars::pid(p) == pvars::kPhoton && p.ancestor_pdg_code == 111);
      if (is_pi0_photon)
      {
        double ke = pvars::ke(p);
        if (ke > best_ke) { best_ke = ke; best_idx = i; }
      }
    }
    return best_idx;
  }
  REGISTER_SELECTOR(leading_pi0_photon, leading_pi0_photon);

  /**
   * @brief Index of the subleading (second-highest-KE) SPINE truth photon
   * with ancestor_pdg_code == 111.  Returns kNoMatch if fewer than 2 found
   * — this is the "missing shower" case.
   **/
  template <class T>
  size_t subleading_pi0_photon(const T& obj)
  {
    size_t lead_idx = leading_pi0_photon(obj);
    size_t best_idx = kNoMatch;
    double best_ke  = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < obj.particles.size(); ++i)
    {
      if (i == lead_idx) continue;
      const auto& p = obj.particles[i];
      bool is_pi0_photon = false;
      if constexpr (std::is_same_v<T, TType>)
        is_pi0_photon = (pvars::pid(p) == pvars::kPhoton && p.ancestor_pdg_code == 111);
      if (is_pi0_photon)
      {
        double ke = pvars::ke(p);
        if (ke > best_ke) { best_ke = ke; best_idx = i; }
      }
    }
    return best_idx;
  }
  REGISTER_SELECTOR(subleading_pi0_photon, subleading_pi0_photon);

} // end namespace selectors::gOre


// ============================================================
//  TrueParticle scope — per-photon properties (SRParticleTruthDLP)
//  Use with selectors: selector = "leading_pi0_photon" or
//                                 "subleading_pi0_photon"
// ============================================================
namespace pvars::gOre
{
  /** @brief KE of the pi0 decay photon [MeV].  pvars::ke uses energy_init. **/
  template <class T>
  double pi0_photon_ke(const T& p)
  {
    return pvars::ke(p);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_ke, pi0_photon_ke);

  /** @brief Initial total energy of the pi0 decay photon [MeV]. **/
  template <class T>
  double pi0_photon_energy_init(const T& p)
  {
    return static_cast<double>(p.energy_init);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_energy_init, pi0_photon_energy_init);

  /** @brief Total energy deposited by the photon in the detector [MeV].
   *  Equivalent to visE in scan_pi0.C (summed over all wire planes). **/
  template <class T>
  double pi0_photon_energy_deposit(const T& p)
  {
    return static_cast<double>(p.energy_deposit);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_energy_deposit, pi0_photon_energy_deposit);

  // ---- Generation / creation position (p.position = tp_gen in scan_pi0.C) ----
  // This is where the photon was created (= the pi0 decay vertex if inside TPC).

  /** @brief Photon creation position x [cm]  (= pi0 decay vertex x if gen inside TPC) **/
  template <class T>
  double pi0_photon_gen_x(const T& p)
  {
    return static_cast<double>(p.position[0]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_gen_x, pi0_photon_gen_x);

  /** @brief Photon creation position y [cm] **/
  template <class T>
  double pi0_photon_gen_y(const T& p)
  {
    return static_cast<double>(p.position[1]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_gen_y, pi0_photon_gen_y);

  /** @brief Photon creation position z [cm] **/
  template <class T>
  double pi0_photon_gen_z(const T& p)
  {
    return static_cast<double>(p.position[2]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_gen_z, pi0_photon_gen_z);

  // ---- TPC start position (p.start_point = tp_start in scan_pi0.C) ----
  // First point in the active TPC volume.  May differ from gen when the photon
  // is created outside the TPC and propagates in.

  /** @brief Photon TPC-entry position x [cm] **/
  template <class T>
  double pi0_photon_start_x(const T& p)
  {
    return static_cast<double>(p.start_point[0]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_start_x, pi0_photon_start_x);

  /** @brief Photon TPC-entry position y [cm] **/
  template <class T>
  double pi0_photon_start_y(const T& p)
  {
    return static_cast<double>(p.start_point[1]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_start_y, pi0_photon_start_y);

  /** @brief Photon TPC-entry position z [cm] **/
  template <class T>
  double pi0_photon_start_z(const T& p)
  {
    return static_cast<double>(p.start_point[2]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_start_z, pi0_photon_start_z);

  // ---- End position (p.end_point = tp_end in scan_pi0.C) ----

  /** @brief Photon end position x [cm] **/
  template <class T>
  double pi0_photon_end_x(const T& p)
  {
    return static_cast<double>(p.end_point[0]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_end_x, pi0_photon_end_x);

  /** @brief Photon end position y [cm] **/
  template <class T>
  double pi0_photon_end_y(const T& p)
  {
    return static_cast<double>(p.end_point[1]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_end_y, pi0_photon_end_y);

  /** @brief Photon end position z [cm] **/
  template <class T>
  double pi0_photon_end_z(const T& p)
  {
    return static_cast<double>(p.end_point[2]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_end_z, pi0_photon_end_z);

  // ---- Direction ----

  /** @brief Photon direction x (unit vector) **/
  template <class T>
  double pi0_photon_dir_x(const T& p)
  {
    return static_cast<double>(p.start_dir[0]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_dir_x, pi0_photon_dir_x);

  /** @brief Photon direction y (unit vector) **/
  template <class T>
  double pi0_photon_dir_y(const T& p)
  {
    return static_cast<double>(p.start_dir[1]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_dir_y, pi0_photon_dir_y);

  /** @brief Photon direction z (unit vector) **/
  template <class T>
  double pi0_photon_dir_z(const T& p)
  {
    return static_cast<double>(p.start_dir[2]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_dir_z, pi0_photon_dir_z);

  // ---- Momentum [MeV/c] ----

  /** @brief Photon momentum x [MeV/c] **/
  template <class T>
  double pi0_photon_px(const T& p)
  {
    return static_cast<double>(p.momentum[0]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_px, pi0_photon_px);

  /** @brief Photon momentum y [MeV/c] **/
  template <class T>
  double pi0_photon_py(const T& p)
  {
    return static_cast<double>(p.momentum[1]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_py, pi0_photon_py);

  /** @brief Photon momentum z [MeV/c] **/
  template <class T>
  double pi0_photon_pz(const T& p)
  {
    return static_cast<double>(p.momentum[2]);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_pz, pi0_photon_pz);

  // ---- Geometry / topology ----

  /** @brief Is the photon fully contained within the active TPC?
   *  1 = contained, 0 = exits the detector. **/
  template <class T>
  double pi0_photon_is_contained(const T& p)
  {
    return static_cast<double>(p.is_contained);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_is_contained, pi0_photon_is_contained);

  /** @brief Trajectory length of the photon in the active TPC [cm].
   *  Zero for photons that do not interact before leaving the detector. **/
  template <class T>
  double pi0_photon_length(const T& p)
  {
    return static_cast<double>(p.length);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_length, pi0_photon_length);

  // ---- SPINE / ancestry flags ----

  /** @brief Is the photon classified as a primary particle by SPINE?
   *  1 = primary, 0 = secondary (e.g. from EM shower). **/
  template <class T>
  double pi0_photon_is_primary(const T& p)
  {
    return static_cast<double>(pvars::primary_classification(p));
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_is_primary, pi0_photon_is_primary);

  /** @brief PDG code of the photon ancestor (should be 111 = pi0 by construction). **/
  template <class T>
  double pi0_photon_ancestor_pdg(const T& p)
  {
    return static_cast<double>(p.ancestor_pdg_code);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_ancestor_pdg, pi0_photon_ancestor_pdg);

  /** @brief G4 track ID of this photon (for debugging / cross-referencing). **/
  template <class T>
  double pi0_photon_track_id(const T& p)
  {
    return static_cast<double>(p.track_id);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_track_id, pi0_photon_track_id);

  /** @brief G4 track ID of the photon's parent (should be the pi0 track ID). **/
  template <class T>
  double pi0_photon_parent_id(const T& p)
  {
    return static_cast<double>(p.parent_id);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::TrueParticle, pi0_photon_parent_id, pi0_photon_parent_id);

} // end namespace pvars::gOre


// ============================================================
//  MCTruth scope — pi0 GENIE primary (SRTrueParticle fields)
//  Full G4 info: gen/start/end positions, wall codes, G4 process
//  codes, generator code, vis_e, nhit, n_daughters.
// ============================================================
namespace vars::gOre
{
  // ---------------------------------------------------------------------------
  // Internal helper: index of the first pi0 (pdg == 111) in obj.prim.
  // Returns std::numeric_limits<size_t>::max() when not found.
  // ---------------------------------------------------------------------------
  template <class T>
  size_t pion_prim_idx(const T& obj)
  {
    for (size_t i = 0; i < obj.prim.size(); ++i)
      if (obj.prim[i].pdg == 111) return i;
    return std::numeric_limits<size_t>::max();
  }

  // ---- Count ----

  /** @brief Number of pi0 photons (SPINE truth, ancestor_pdg_code == 111)
   *  in this interaction.  Expects 2 for a normal pi0 → γγ decay. **/
  template <class T>
  double n_pi0_photons(const T& obj)
  {
    size_t n = 0;
    for (const auto& p : obj.particles)
      if (pvars::pid(p) == pvars::kPhoton && p.ancestor_pdg_code == 111) ++n;
    return static_cast<double>(n);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::True, n_pi0_photons, n_pi0_photons);

  // ---- MCTruth: Number of pi0s in GENIE primary list ----

  /** @brief Number of pi0s in the GENIE primary list. **/
  template <class T>
  double n_genie_pions(const T& obj)
  {
    size_t n = 0;
    for (size_t i = 0; i < obj.prim.size(); ++i)
      if (obj.prim[i].pdg == 111) ++n;
    return static_cast<double>(n);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, n_genie_pions, n_genie_pions);

  /** @brief Pi0 total energy at generation [MeV]  (genE × 1000). **/
  template <class T>
  double pion_gen_energy(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].genE) * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_energy, pion_gen_energy);

  /** @brief Pi0 total energy at its first point in the active TPC [MeV]. **/
  template <class T>
  double pion_start_energy(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].startE) * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_start_energy, pion_start_energy);

  // ---- Generation position ----

  /** @brief Pi0 generation position x [cm] **/
  template <class T>
  double pion_gen_x(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].gen.x);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_x, pion_gen_x);

  /** @brief Pi0 generation position y [cm] **/
  template <class T>
  double pion_gen_y(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].gen.y);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_y, pion_gen_y);

  /** @brief Pi0 generation position z [cm] **/
  template <class T>
  double pion_gen_z(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].gen.z);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_z, pion_gen_z);

  // ---- Generation momentum [MeV/c] ----

  /** @brief Pi0 momentum x at generation [MeV/c] **/
  template <class T>
  double pion_gen_px(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].genp.x) * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_px, pion_gen_px);

  /** @brief Pi0 momentum y at generation [MeV/c] **/
  template <class T>
  double pion_gen_py(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].genp.y) * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_py, pion_gen_py);

  /** @brief Pi0 momentum z at generation [MeV/c] **/
  template <class T>
  double pion_gen_pz(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].genp.z) * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_gen_pz, pion_gen_pz);

  // ---- Start position (first point in active TPC) ----

  /** @brief Pi0 start position x in the active TPC volume [cm] **/
  template <class T>
  double pion_start_x(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].start.x);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_start_x, pion_start_x);

  /** @brief Pi0 start position y in the active TPC volume [cm] **/
  template <class T>
  double pion_start_y(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].start.y);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_start_y, pion_start_y);

  /** @brief Pi0 start position z in the active TPC volume [cm] **/
  template <class T>
  double pion_start_z(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].start.z);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_start_z, pion_start_z);

  // ---- End position (pi0 decay vertex) ----

  /** @brief Pi0 end / decay vertex x [cm] **/
  template <class T>
  double pion_end_x(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].end.x);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_end_x, pion_end_x);

  /** @brief Pi0 end / decay vertex y [cm] **/
  template <class T>
  double pion_end_y(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].end.y);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_end_y, pion_end_y);

  /** @brief Pi0 end / decay vertex z [cm] **/
  template <class T>
  double pion_end_z(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].end.z);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_end_z, pion_end_z);

  // ---- Geometry flags ----

  /** @brief Is the pi0 fully contained within the active TPC?  1=yes 0=no. **/
  template <class T>
  double pion_contained(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].contained);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_contained, pion_contained);

  /** @brief Wall the pi0 exits through (Wall_t enum, see file header). **/
  template <class T>
  double pion_wallout(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].wallout);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_wallout, pion_wallout);

  /** @brief Wall the pi0 enters through (Wall_t enum).
   *  0 = kWallNone (created inside active volume). **/
  template <class T>
  double pion_wallin(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].wallin);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_wallin, pion_wallin);

  // ---- G4 process codes ----

  /** @brief G4 process at pi0 creation (g4_process_ enum).
   *  Expected: 0 (kG4primary) — created by GENIE at the neutrino vertex. **/
  template <class T>
  double pion_start_process(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].start_process);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_start_process, pion_start_process);

  /** @brief G4 process at pi0 end (g4_process_ enum).
   *  Expected: 3 (kG4Decay) — pi0 → γγ. **/
  template <class T>
  double pion_end_process(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].end_process);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_end_process, pion_end_process);

  /** @brief Generator that created the pi0 (generator_ enum).
   *  Expected: 1 (kGENIE). 0 = kUnknownGenerator (empty). **/
  template <class T>
  double pion_generator(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].generator);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_generator, pion_generator);

  // ---- Trajectory length and deposited energy ----

  /** @brief Pi0 trajectory length in the active TPC [cm]. **/
  template <class T>
  double pion_length(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].length);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_length, pion_length);

  /** @brief Total visible energy from the pi0 across all planes [MeV].
   *  Expected 0 for a neutral pi0 (Dalitz decay would be non-zero). **/
  template <class T>
  double pion_vis_e(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    const auto& p = obj.prim[idx];
    double tot = 0.;
    for (int cryo = 0; cryo < 2; ++cryo)
      for (int pl = 0; pl < 3; ++pl)
        tot += static_cast<double>(p.plane[cryo][pl].visE);
    return tot * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_vis_e, pion_vis_e);

  /** @brief Total number of wire hits from the pi0 across all planes.
   *  Expected 0 for a neutral pi0. **/
  template <class T>
  double pion_nhit(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    const auto& p = obj.prim[idx];
    unsigned tot = 0;
    for (int cryo = 0; cryo < 2; ++cryo)
      for (int pl = 0; pl < 3; ++pl)
        tot += p.plane[cryo][pl].nhit;
    return static_cast<double>(tot);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_nhit, pion_nhit);

  /** @brief Number of G4 daughters of the pi0.
   *  pi0 → γγ gives 2; Dalitz pi0 → γe+e− gives 3. **/
  template <class T>
  double pion_n_daughters(const T& obj)
  {
    size_t idx = pion_prim_idx(obj);
    if (idx == std::numeric_limits<size_t>::max())
      return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(obj.prim[idx].daughters.size());
  }
  REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_n_daughters, pion_n_daughters);

} // end namespace vars::gOre


// ============================================================
//  G4Truth scope — pi0 decay photon properties from rec.true_particles
//  Mirrors rec.true_particles branches from scan_pi0.C lines 208-245.
//  Variables are named photon_leading_* / photon_subleading_* where
//  "leading" = higher startE, "subleading" = lower startE.
//  TOML usage:  { name = "photon_leading_g4id", type = "g4truth" }
//  Branch names: g4truth_photon_leading_g4id, g4truth_photon_subleading_g4id, ...
//
//  G4TruthContext:  ctx.mc = SRTrueInteraction (GENIE truth)
//                   ctx.sr = StandardRecord    (for true_particles flat array)
// ============================================================
namespace vars::gOre
{
  namespace g4detail
  {
    /// G4ID of the first pi0 (pdg==111) in ctx.mc->prim. Returns -1 if none.
    inline int pi0_g4id(const G4TruthType& ctx)
    {
      if (!ctx.mc) return -1;
      for (size_t i = 0; i < ctx.mc->prim.size(); ++i)
        if (ctx.mc->prim[i].pdg == 111)
          return static_cast<int>(ctx.mc->prim[i].G4ID);
      return -1;
    }

    /// Index into ctx.sr->true_particles of the rank-th highest-startE photon
    /// whose parent is the pi0.  rank 0 = leading, rank 1 = subleading.
    /// Returns -1 if not found.
    inline int64_t photon_idx(const G4TruthType& ctx, int rank)
    {
      int pi0id = pi0_g4id(ctx);
      if (pi0id < 0 || !ctx.sr) return -1;
      std::vector<std::pair<double, int64_t>> photons;
      for (int64_t i = 0; i < static_cast<int64_t>(ctx.sr->true_particles.size()); ++i)
      {
        const auto& tp = ctx.sr->true_particles[i];
        if (static_cast<int>(tp.parent) == pi0id && tp.pdg == 22)
          photons.emplace_back(static_cast<double>(tp.startE), i);
      }
      std::sort(photons.begin(), photons.end(),
                [](const auto& a, const auto& b){ return a.first > b.first; });
      if (rank < static_cast<int>(photons.size())) return photons[rank].second;
      return -1;
    }
  } // namespace g4detail

// Macro: generates leading (rank 0) and subleading (rank 1) variable pairs.
// EXPR may reference `tp` (const auto& = ctx.sr->true_particles[idx]).
#define G4_PHOTON_VAR_IMPL(RANKNAME, RANKNUM, FNAME, EXPR)               \
  template <class T>                                                      \
  double photon_##RANKNAME##_##FNAME(const T& ctx)                       \
  {                                                                       \
    int64_t idx = g4detail::photon_idx(ctx, RANKNUM);                    \
    if (idx < 0) return kNoMatchValue;                                    \
    const auto& tp = ctx.sr->true_particles[idx];                        \
    return (EXPR);                                                        \
  }                                                                       \
  REGISTER_VAR_SCOPE(RegistrationScope::G4Truth,                         \
      photon_##RANKNAME##_##FNAME, photon_##RANKNAME##_##FNAME);

#define G4_PHOTON_VAR(FNAME, EXPR)                 \
  G4_PHOTON_VAR_IMPL(leading,    0, FNAME, EXPR)   \
  G4_PHOTON_VAR_IMPL(subleading, 1, FNAME, EXPR)

  // ---- Identification ----
  G4_PHOTON_VAR(g4id,           static_cast<double>(tp.G4ID))
  G4_PHOTON_VAR(parent,         static_cast<double>(tp.parent))
  G4_PHOTON_VAR(pdg,            static_cast<double>(tp.pdg))
  G4_PHOTON_VAR(interaction_id, static_cast<double>(tp.interaction_id))

  // ---- Energy [MeV] — startE stored in GeV, ×1000 ----
  G4_PHOTON_VAR(start_e,        static_cast<double>(tp.startE) * 1000.)

  // ---- Generation (creation) position [cm] ----
  G4_PHOTON_VAR(gen_x,          static_cast<double>(tp.gen.x))
  G4_PHOTON_VAR(gen_y,          static_cast<double>(tp.gen.y))
  G4_PHOTON_VAR(gen_z,          static_cast<double>(tp.gen.z))

  // ---- TPC start position [cm] ----
  G4_PHOTON_VAR(start_x,        static_cast<double>(tp.start.x))
  G4_PHOTON_VAR(start_y,        static_cast<double>(tp.start.y))
  G4_PHOTON_VAR(start_z,        static_cast<double>(tp.start.z))

  // ---- End position [cm] ----
  G4_PHOTON_VAR(end_x,          static_cast<double>(tp.end.x))
  G4_PHOTON_VAR(end_y,          static_cast<double>(tp.end.y))
  G4_PHOTON_VAR(end_z,          static_cast<double>(tp.end.z))

  // ---- Start momentum [MeV/c] — startp stored in GeV/c, ×1000 ----
  G4_PHOTON_VAR(startp_x,       static_cast<double>(tp.startp.x) * 1000.)
  G4_PHOTON_VAR(startp_y,       static_cast<double>(tp.startp.y) * 1000.)
  G4_PHOTON_VAR(startp_z,       static_cast<double>(tp.startp.z) * 1000.)

  // ---- Geometry ----
  G4_PHOTON_VAR(length,         static_cast<double>(tp.length))
  G4_PHOTON_VAR(contained,      static_cast<double>(tp.contained))
  G4_PHOTON_VAR(wallout,        static_cast<double>(tp.wallout))
  G4_PHOTON_VAR(wallin,         static_cast<double>(tp.wallin))

  // ---- G4 process / generator codes ----
  G4_PHOTON_VAR(start_process,  static_cast<double>(tp.start_process))
  G4_PHOTON_VAR(end_process,    static_cast<double>(tp.end_process))
  G4_PHOTON_VAR(generator,      static_cast<double>(tp.generator))

#undef G4_PHOTON_VAR_IMPL
#undef G4_PHOTON_VAR

  // ---- Deposited energy [MeV] and wire hits — summed over all 6 planes ----
  // plane[cryo][wire_plane], cryo in {0,1}, wire_plane in {0,1,2}
  // visE stored in GeV, ×1000 converted to MeV.

  template <class T>
  double photon_leading_vis_e(const T& ctx)
  {
    int64_t idx = g4detail::photon_idx(ctx, 0);
    if (idx < 0) return kNoMatchValue;
    const auto& tp = ctx.sr->true_particles[idx];
    double tot = 0.;
    for (int c = 0; c < 2; ++c)
      for (int p = 0; p < 3; ++p)
        tot += static_cast<double>(tp.plane[c][p].visE);
    return tot * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::G4Truth, photon_leading_vis_e, photon_leading_vis_e);

  template <class T>
  double photon_subleading_vis_e(const T& ctx)
  {
    int64_t idx = g4detail::photon_idx(ctx, 1);
    if (idx < 0) return kNoMatchValue;
    const auto& tp = ctx.sr->true_particles[idx];
    double tot = 0.;
    for (int c = 0; c < 2; ++c)
      for (int p = 0; p < 3; ++p)
        tot += static_cast<double>(tp.plane[c][p].visE);
    return tot * 1000.;
  }
  REGISTER_VAR_SCOPE(RegistrationScope::G4Truth, photon_subleading_vis_e, photon_subleading_vis_e);

  template <class T>
  double photon_leading_nhit(const T& ctx)
  {
    int64_t idx = g4detail::photon_idx(ctx, 0);
    if (idx < 0) return kNoMatchValue;
    const auto& tp = ctx.sr->true_particles[idx];
    unsigned tot = 0;
    for (int c = 0; c < 2; ++c)
      for (int p = 0; p < 3; ++p)
        tot += tp.plane[c][p].nhit;
    return static_cast<double>(tot);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::G4Truth, photon_leading_nhit, photon_leading_nhit);

  template <class T>
  double photon_subleading_nhit(const T& ctx)
  {
    int64_t idx = g4detail::photon_idx(ctx, 1);
    if (idx < 0) return kNoMatchValue;
    const auto& tp = ctx.sr->true_particles[idx];
    unsigned tot = 0;
    for (int c = 0; c < 2; ++c)
      for (int p = 0; p < 3; ++p)
        tot += tp.plane[c][p].nhit;
    return static_cast<double>(tot);
  }
  REGISTER_VAR_SCOPE(RegistrationScope::G4Truth, photon_subleading_nhit, photon_subleading_nhit);

} // end namespace vars::gOre (G4Truth section)

#endif // PION_VARS_GORE_H
