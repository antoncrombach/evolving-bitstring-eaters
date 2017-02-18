//
// Part of the observer/subject pattern
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "subject.hh"
#include "observer.hh"
#include "stream_manager.hh"

void 
dorsalfin::Subject::attach( Observer *o ) {
    obs_.push_back( o );
}

void 
dorsalfin::Subject::detach( Observer *o ) {
    // observer 'o' exists in obs_
    obs_.erase( std::find( obs_.begin(), obs_.end(), o ) );
}

void 
dorsalfin::Subject::detachAll() {
    obs_.clear();
}

void 
dorsalfin::Subject::notify() 
{ notify( this ); }

void
dorsalfin::Subject::notify( Subject *s ) {
    for( obs_iter i = obs_.begin(); i != obs_.end(); ++i ) {
        ( **i ).update( s );
    }
}

