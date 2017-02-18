//
// Builds factory, population, notifies observers.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_MODEL_H_
#define _DORSALFIN_MODEL_H_

#include "defs.hh"
#include "fluke.hh"
#include "factory.hh"
#include "event.hh"

namespace dorsalfin {

    /// \class Model
    /// \brief Build the model, then let the population and environment run
    /// for some time.
    ///
    /// The model consists of a population of agents and an environment. Both
    /// are assembled in the factory according to the given configuration. 
    /// The observers keep track of certain features of the agents and the 
    /// state of the environment.
    class Model {
        public:
            /// Constructor
            Model( Fluke * );
            /// Destructor
            ~Model();

            /// Build the model. This includes building a population of
            /// agents, creating an observer manager and constructing an
            /// environment.
            void build();
            /// Initialise all the components
            void initialise();
            /// Rebuild / configure for a new run
            void rebuild();
            /// Perform one timestep
            void step();
            /// Round-up of the simulation
            void finish();

            /// What time is it?
            long now();
            /// Did the simulation end yet?
            bool hasEnded();
            
            /// Return the observer manager
            ObserverManager & observerManager() const;
            /// Return the population
            Population & population() const;
            /// Return the environment
            Environment & environment() const;

        public:
            /// Set the ending time of the simulation run
            static void endTime( long );
            /// Get the ending time.
            static long endTime();

        private:
            /// Build observers and connect them to the right subject
            void observe();
            /// Destroy observers and close their output streams.
            void unobserve();
            
        private:
            static long end_time_;

        private:
            long time_;
            
            Factory factory_;
            Fluke *fluke_;
            Population *poppy_, *cache_poppy_;
            Environment *environ_;
            ObserverManager *observers_;
            
            std::priority_queue< Event* > *events_;
    };

    inline ObserverManager & Model::observerManager() const
    { return *observers_; }

    inline Population & Model::population() const
    { return *poppy_; }
    
    inline Environment & Model::environment() const
    { return *environ_; }
    
    inline long Model::now()
    { return time_; }

    inline void Model::endTime( long t )
    { end_time_ = t; }

    inline long Model::endTime()
    { return end_time_; }
}
#endif

