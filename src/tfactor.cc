//
// Transcription factor, as special type of gene
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "tfactor.hh"
#include "pool.hh"

template<> dorsalfin::ObjectCache< dorsalfin::TransFactor >* 
dorsalfin::ObjectCache< dorsalfin::TransFactor >::instance_ = 0;


dorsalfin::TransFactor::TransFactor() : Gene() {}
dorsalfin::TransFactor::TransFactor( int t ) : Gene( t ) {}
dorsalfin::TransFactor::TransFactor( int t, int s, int th ) : Gene( t, s, th ) {}
dorsalfin::TransFactor::TransFactor( const TransFactor &ds ) :  Gene( ds ) {}

dorsalfin::ChromosomeElement* 
dorsalfin::TransFactor::clone() const {
    //return new TransFactor( *this );
    TransFactor *g = 
        ObjectCache< TransFactor >::instance()->borrowObject();
    g->copy( *this );
    return g;
}

bool
dorsalfin::TransFactor::toPool() {
    return ObjectCache< TransFactor >::instance()->returnObject( this );
}

// one method overriden
void 
dorsalfin::TransFactor::writeXml( std::ostream &os ) const { 
    os << "<tfac id=\"" + boost::lexical_cast< std::string >( tag_ ) +
    "\" st=\"" + boost::lexical_cast< std::string >( state_ ) +
    "\" th=\"" + boost::lexical_cast< std::string >( threshold_ ) + "\"/>\n";
}

