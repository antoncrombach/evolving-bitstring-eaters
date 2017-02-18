//
// Manages observers that collect stats and log them to file.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_COLLECTING_H_
#define _DORSALFIN_COLLECTING_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class ObserverManager
    /// \brief Mediator connecting subjects to observers.
    ///
    /// In order to decouple the binding between observers and their subjects,
    /// and to encapsulate the responsibilities of observers, we introduce a
    /// mediator class that manages the subscription and notification of
    /// both observers and subjects.
    class ObserverManager {
        public:
            /// Iterator for observers that log periodically to a stream
            typedef std::multimap< Subject*, LogObserver* >::iterator
                sub_obs_iter;
            /// Iterator for observsers that log if necessary (when changes
            /// occur)
            typedef std::multimap< Subject*, AsyncLogObserver* >::iterator
                sub_asyn_obs_iter;

        public:
            /// Constructor
            ObserverManager();
            /// Destructor
            ~ObserverManager();

            /// Subscribe a periodic observer to a subject
            void subscribe( Subject*, LogObserver* );
            /// Subscribe an asynchronous observer to a subject
            void subscribe( Subject*, AsyncLogObserver* );
            /// Remove the subscription of an observer from its subject
            void unsubscribe( Subject*, LogObserver* );
            /// Remove all subscriptions of this subject
            void unsubscribe( Subject* );
            /// Remove all subscription of the given observer
            void unsubscribe( LogObserver* );

            /// Notify all observers of this subject (only periodic observers
            /// are affected)
            void notify( Subject* );
            /// Notify all periodic observers
            void notifyAll(); 
            /// Close the logs of all observers
            void closeAll();

        private:
            std::multimap< Subject*, LogObserver* > subobs_;
            std::multimap< Subject*, AsyncLogObserver* > subasynobs_;
    };
}
#endif

