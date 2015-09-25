# SaRoMaN
Simulation and reconstruction of muons and neutrinos

Comprehensive software for MIND/nuSTORM simulations 


The intention is to re-evaluate and rededicate the existing simulation
software associated with nuSTORM using the lastest versions of GENIE,
GEANT4, and RecPack. The motivation for this is that the current
version of the simulation software is either heavily cludged (in the
case of the MIND simulation) or it does not include a fully functional
reconstruction (as is the case of gnomon and potentially LBNOsoft).

The goal is to make a modular software package similar to what
currently exists but with some expanded functionality:

* use of geometry in GENIE allowing for the simulation of detectors by
  neutrino rate (rather than interaction rate) and the simulation of
  near detectors for the purpose of near-far extrapolations

* baseline integration of CRY cosmic ray simulation.

* default use of GDML geometry definition for increased flexibility in
  geometry. This should make the scintillator segmentation of
  SuperBIND trivial.

* extraction of GDML geometry in reconstruction application for
  simplified RecPack model.

* sensible digitization that works with the flexible geometry. This
  should reduce sensible groupings of ionization into signals that
  mimic SiPM response. Clustering should be moved to reconstruction.

* python binding/control of entire (or partial) simulation cycle to
  streamline the user experiance.

The prime concern is the use of a sensible data format through the
entire framework --- the only problem with assuming Yordan's LBNO
framework whole is getting it to communicate with RecPack. The
feasibility of the communication of hits in this format needs to be
evaluated.

List of external packages required;
* GENIE
* GEANT4
* RecPack
* bhep?
* LHAPDF
* Xerces
* expat
* others?

The workspace should be composed of the following packages

* GENDAfS -- Genie Enabled Neutrino Detector Analogue for
  Simulation. -- This is where the simulation is actually
  implemented. Events produced with GENIE are transported through the
  GEANT4 generated detector. The energy deposited in each detector
  element is recorded and passed to digitization

* DoSIiS -- Digitization of Simulated Intractions in Scintillator. --
  Energy deposition in each detector element is attenuated and smeared
  as in a real device. Modelling of response of detector elements such
  as scintillator bars and SiPMs should exist here. Output should be
  number of photo-electrons collected per detector element.

* RoNDaGoST -- Reconstruction of Neutrino Detector and Generation of
  Software Trees. -- The digitized output is rendered into hit
  positions that can be used in particle identification and track
  reconstruction. Particle identification and fitting based on the
  Kalman filters and fitters provided by RecPack.

* AoRMaS -- Analysis of Reconstructed Muons and Showers. -- A bundle
  of methods used to select muons, or isolate showers, from
  reconstructed events. Should have a cuts based option as well as an
  option that uses multivariate methods. Can be streamlined from
  currently existing software.