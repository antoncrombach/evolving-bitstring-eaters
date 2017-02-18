//
// Observer/subject pattern (GOF).
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_OBSERVER_H_
#define _DORSALFIN_OBSERVER_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class Observer
    /// \brief Abstract base class of observer/subject pattern (GOF).
    ///
    /// The observer/subject pattern is used to decouple the actual model from
    /// data gathering (separation of concerns). In \c fluke the observers 
    /// collect data and statistics from the model.
    class Observer {
        public:
        /// Constructor
        Observer() {};
        /// Copy constructor
        Observer( const Observer &o ) {};
        /// Destructor
        virtual ~Observer() {};

        /// Starting the logging, any initialization needed?
        virtual void initialize( Subject* ) = 0;
        /// Update \c this state according to the changes in \c Subject
        virtual void update( Subject* ) = 0;
        /// Ending the logging, like writing a footer.
        virtual void finalize() = 0;        
    };

    /// \class LogObserver
    /// \brief Observer specialised in logging collected data in a file.
    class LogObserver : public Observer {
        public:
        /// Constructor. The integer specifies the period of logging 
        /// (in simulation timesteps).
        LogObserver( long );
        /// Copy constructor
        LogObserver( const LogObserver & );
        /// Destructor
        virtual ~LogObserver();

        /// Open a log (file) stream
        void openLog( std::string );
        /// Close the log stream
        void closeLog();

        /// Update with intervals of a given period. We applied 
        /// a template-method pattern (GOF) here.
        virtual void update( Subject* );
        /// Actual updating. This function needs to be overriden in child
        /// classes.
        virtual void doUpdate( Subject* ) = 0;

        protected:
        /// Log stream to write to
        boost::filesystem::ofstream *log_;
        /// Period of writing (in simulation timesteps)
        long interval_, val_;
    };

    /// \class AsyncLogObserver
    /// \brief Observer specialised in asynchronously logging to file
    class AsyncLogObserver : public Observer {
        public:
        /// Constructor.
        AsyncLogObserver();
        /// Copy constructor
        AsyncLogObserver( const AsyncLogObserver & );
        /// Destructor
        virtual ~AsyncLogObserver();

        // update needs to be overloaded by child classes
        /// Open a log (file) stream
        void openLog( std::string );
        /// Close the log stream
        void closeLog();

        protected:
        /// Log stream to write to
        boost::filesystem::ofstream *log_;
    };
}
#endif

