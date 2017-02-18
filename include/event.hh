//
// Simple system for handling discrete events
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_EVENT_H_
#define _DORSALFIN_EVENT_H_

#include "defs.hh"

namespace dorsalfin {
    /// \class Event
    class Event {
        public:
        /// dtor
        virtual ~Event() {}

        /// What happens...
        virtual void execute( Population & ) = 0;
        /// When is this event supposed to happen?
        long when() const;

        protected:
        /// ctor
        Event() : time_( 0 ) {}
        /// ctor, set time
        Event( long t ) : time_( t ) {}
        
        protected:
        long time_;
    };
    
    inline long Event::when() const
    { return time_; }
    
    inline bool operator<( const Event &x, const Event &y )
    { return x.when() > y.when(); }
    
    /// \class ZeroRatesEvent
    class ZeroRatesEvent : public Event {
        public:
        /// ctor
        ZeroRatesEvent( long t ) : Event( t ) {}
        /// dtor
        virtual ~ZeroRatesEvent() {}
        
        virtual void execute( Population & );
    };
}
#endif
