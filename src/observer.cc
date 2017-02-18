//
// Implementation of a few methods in the subject/observer pattern.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "observer.hh"
#include "stream_manager.hh"


dorsalfin::LogObserver::LogObserver( long i )
    : Observer(), interval_( i ), val_( 1 ) {
    log_ = 0;
}

dorsalfin::LogObserver::LogObserver( const LogObserver &lo ) 
    : Observer( lo ) {
    log_ = lo.log_;
    interval_ = lo.interval_;
    val_ = lo.val_;
}

dorsalfin::LogObserver::~LogObserver() {
    closeLog();
}

void 
dorsalfin::LogObserver::openLog( std::string fname ) {
    log_ = StreamManager::instance()->openOutFileStream( fname, 
            std::fstream::out | std::fstream::app );
}

void 
dorsalfin::LogObserver::closeLog() {
    if( log_ != 0 ) {
        finalize();
        StreamManager::instance()->closeOutFileStream( log_ );
        log_ = 0;
    }
}

void
dorsalfin::LogObserver::update( Subject *s ) {
    // template method
    if( val_ == 1 ) {
        val_ = interval_;
        doUpdate( s );
    } else {
        --val_;
    }
}

//
// async log observer
//
dorsalfin::AsyncLogObserver::AsyncLogObserver() : Observer() {
    log_ = 0;
}

dorsalfin::AsyncLogObserver::AsyncLogObserver( const AsyncLogObserver &lo ) 
    : Observer( lo ) {
    log_ = lo.log_;
}

dorsalfin::AsyncLogObserver::~AsyncLogObserver() {
    closeLog();
}

void 
dorsalfin::AsyncLogObserver::openLog( std::string fname ) {
    log_ = StreamManager::instance()->openOutFileStream( fname, 
            std::fstream::out | std::fstream::app );
}

void 
dorsalfin::AsyncLogObserver::closeLog() {
    if( log_ != 0 ) {
        finalize();
        StreamManager::instance()->closeOutFileStream( log_ );
        log_ = 0;
    }
}

