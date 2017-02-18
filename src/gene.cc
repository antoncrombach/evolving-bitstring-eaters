//
// Gene with interactions and threshold?
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "gene.hh"
#include "pool.hh"

template<> dorsalfin::ObjectCache< dorsalfin::Gene >* 
dorsalfin::ObjectCache< dorsalfin::Gene >::instance_ = 0;


dorsalfin::Gene::Gene() 
    : ChromosomeElement(), tag_( 0 ), state_( 0 ), 
      threshold_( 0 ) {}

dorsalfin::Gene::Gene( int t ) 
    : ChromosomeElement(), tag_( t ), state_( 0 ),
      threshold_( 0 ) {}

dorsalfin::Gene::Gene( int t, int s, int th ) 
    : ChromosomeElement(), tag_( t ), state_( s ), threshold_( th ) {}

dorsalfin::Gene::Gene( const Gene &ds ) 
    :  ChromosomeElement( ds ), tag_( ds.tag_ ), state_( ds.state_ ),
       threshold_( ds.threshold_ ) {}

dorsalfin::ChromosomeElement* 
dorsalfin::Gene::clone() const {
    //return new Gene( *this );
    Gene *g = ObjectCache< Gene >::instance()->borrowObject();
    g->copy( *this );
    return g;
}
      
void
dorsalfin::Gene::copy( const ChromosomeElement &ce ) {
    const Gene &g = dynamic_cast< const Gene& >( ce );
    tag_ = g.tag_;
    state_ = g.state_;
    threshold_ = g.threshold_;
}

bool
dorsalfin::Gene::toPool() {
    return ObjectCache< Gene >::instance()->returnObject( this );
}

void
dorsalfin::Gene::mutateThreshold() {
    // only five values {-2,-1,0,1,2}
    threshold_ = static_cast< int >( uniform() * 4 + 0.5 ) - 2;
}

void
dorsalfin::Gene::mutateState() {
    state_ = 1 - state_;
}
void
dorsalfin::Gene::import( const RefNode &r ) {
    //tag_ = r.node_id;
    //state_ = 0;
    threshold_ = r.threshold;
}

void 
dorsalfin::Gene::writeXml( std::ostream &os ) const { 
    os << "<gene id=\"" + boost::lexical_cast< std::string >( tag_ ) +
    "\" st=\"" + boost::lexical_cast< std::string >( state_ ) +
    "\" th=\"" + boost::lexical_cast< std::string >( threshold_ ) + "\"/>\n";
}

