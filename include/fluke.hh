//
// The heart of the program...
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_MAIN_H_
#define _DORSALFIN_MAIN_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class Fluke
    /// \brief The heart of the program (or maybe the brain?)
    ///
    /// Fluke is a small class that connects the model to the configuration
    /// and the filestream manager. It takes care of initialisation of a run,
    /// makes sure everything is set up correctly (aborts on errors) and
    /// runs the simulation.
    class Fluke {
        public:
        /// Constructor needing \c argc and \c argv to parse command line
        /// arguments.
        Fluke( int, char** );
        /// Destructor
        ~Fluke();

        /// Run the simulation(s)
        int run();
        
        /// Return the configuration.
        Config & configuration() const;
        /// Return the model.
        Model & model() const;

        private:
        // run the simulations        
        void simulate();
        // run the sim
        void doRun();
        // after a run, reconfigure for next run
        void reconfigure( int );

        private:
        Config *config_;
        Model *model_;
        std::string config_fname_;
    };

    inline Config & Fluke::configuration() const
    { return *config_; }

    inline Model & Fluke::model() const
    { return *model_; }
}
#endif

