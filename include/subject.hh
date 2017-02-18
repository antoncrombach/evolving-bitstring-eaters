//
// !one-liner about file contents
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_SUBJECT_H_
#define _DORSALFIN_SUBJECT_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class Subject
    /// \brief Abstract base class of observer/subject pattern (GOF).
    ///
    /// The observer/subject pattern is used to decouple the actual model from
    /// data gathering (separation of concerns). In \c fluke the subjects are 
    /// parts of the model, f.i. the population, a specific agent, genes.
    class Subject {
        public:
            /// Observer list iterator
            typedef std::list< Observer* >::iterator obs_iter;

        public:
            /// Constructor
            Subject() : obs_() {};
            /// Copy constructor
            Subject( const Subject &s ) : obs_( s.obs_ ) {};
            /// Destructor
            virtual ~Subject() {};

            /// Register an observer with \c this
            void attach( Observer* );
            /// Unregister an observer with \c this
            /// \pre observer has been registered
            void detach( Observer* );
            /// Unregister all observers
            void detachAll();

            /// Notify all observers (because they need to update)
            void notify();
            /// Notify the observers with any subject
            void notify( Subject * );

        private:
            std::list< Observer* > obs_;
    };
}
#endif
