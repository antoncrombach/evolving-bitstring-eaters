//
// Implementation of the observer manager.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "observer_manager.hh"
#include "observer.hh"
#include "subject.hh"

dorsalfin::ObserverManager::ObserverManager() : subobs_(), subasynobs_() {}

dorsalfin::ObserverManager::~ObserverManager() {
    for( sub_obs_iter i = subobs_.begin(); i != subobs_.end(); ++i ) {
        delete i->second;
    }
    for( sub_asyn_obs_iter i = subasynobs_.begin();
        i != subasynobs_.end(); ++i ) {
        delete i->second;
    }
    subobs_.clear();
    subasynobs_.clear();
}

void 
dorsalfin::ObserverManager::subscribe( Subject *s, LogObserver *lo ) {
    lo->initialize( s );
    subobs_.insert( std::make_pair( s, lo ) );
}

void 
dorsalfin::ObserverManager::subscribe( Subject *s, AsyncLogObserver *lo ) {
    s->attach( lo );
    lo->initialize( s );
    subasynobs_.insert( std::make_pair( s, lo ) );
}

void 
dorsalfin::ObserverManager::unsubscribe( Subject *s, LogObserver *lo ) {
    std::pair< sub_obs_iter, sub_obs_iter > aux =
        subobs_.equal_range( s );
    sub_obs_iter i = aux.first;
    while( i != aux.second ) {
        if( i->second == lo ) {
            sub_obs_iter j = i;
            ++i;
            subobs_.erase( j );
        } else {
            ++i;
        }
    }
}

void 
dorsalfin::ObserverManager::unsubscribe( Subject *s ) {
    std::pair< sub_obs_iter, sub_obs_iter > aux =
        subobs_.equal_range( s );
    subobs_.erase( aux.first, aux.second );
}

void 
dorsalfin::ObserverManager::unsubscribe( LogObserver *lo ) {
    // unfortunately we have to iterate through the entire map
    sub_obs_iter i = subobs_.begin(); 
    while( i != subobs_.end() ) {
        if( i->second == lo ) {
            sub_obs_iter j = i;
            ++i;
            subobs_.erase( j );
        } else {
            ++i;
        }
    }
}

void 
dorsalfin::ObserverManager::notify( Subject *s ) {
    std::pair< sub_obs_iter, sub_obs_iter > aux = 
        subobs_.equal_range( s );
    for( sub_obs_iter i = aux.first; i != aux.second; ++i ) {
         i->second->update( i->first );
    }
}

void 
dorsalfin::ObserverManager::notifyAll() {
    for( sub_obs_iter i = subobs_.begin(); i != subobs_.end(); ++i ) {
        i->second->update( i->first );
    }
}

void 
dorsalfin::ObserverManager::closeAll() {
    for( sub_obs_iter i = subobs_.begin(); i != subobs_.end(); ++i ) {
        i->second->closeLog();
    }
    for( sub_asyn_obs_iter i = subasynobs_.begin(); 
        i != subasynobs_.end(); ++i ) {
        i->second->closeLog();
    }
}
