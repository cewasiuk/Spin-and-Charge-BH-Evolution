// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_PARTICLE_COUNTS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_PARTICLE_COUNTS);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      foreach(const long & id , _ids) {
	stringstream ss1;
	ss1 << "Hist" << id;
	_histograms[id] = bookHisto1D(ss1.str(), 100,-10.,0.);
	stringstream ss2;
	ss2 << "Count" << id;
	_counts[id] = bookCounter(ss2.str());
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double ebeam = 0.5*sqrtS();
      double weight = event.weight();
      const Particles& particles = applyProjection<FinalState>(event, "FS").particles();
      foreach(Particle p, particles) {
	long id = abs(p.pdgId());
	if(_ids.find(id)==_ids.end()) {
	 /// cerr << "Unknown type of stable particle "
	 ///      << p << "\n";
	  continue;
	}
	_histograms[id]->fill(log(p.momentum().E()/ebeam),weight);
	_counts[id]->fill(weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize(_h_YYYY); // normalize to unity
      // scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      foreach(const long & id , _ids) {
	scale(_histograms[id], 1./sumOfWeights());
	scale(_counts[id], 1./sumOfWeights());
      }
    }

    //@}


    /// @name Histograms
    //@{
    map<long,Histo1DPtr> _histograms;
    map<long,CounterPtr> _counts;
    set<long> _ids={211,111,22,12,14,16,130,321,2212,2112,11,13};
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_PARTICLE_COUNTS);


}
