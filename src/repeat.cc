//
// Implementation of a special chromosome element, the repeat.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "repeat.hh"
#include "pool.hh"

template<> dorsalfin::ObjectCache< dorsalfin::Repeat >* 
dorsalfin::ObjectCache< dorsalfin::Repeat >::instance_ = 0;

dorsalfin::Repeat::Repeat() : ChromosomeElement(), dsb_( false ) {}

dorsalfin::Repeat::Repeat( const Repeat &tp ) : ChromosomeElement( tp ) {
    Repeat::copy( tp );
}

dorsalfin::ChromosomeElement* 
dorsalfin::Repeat::clone() const {
    //return new Repeat( *this );
    Repeat *ltr = ObjectCache< Repeat >::instance()->borrowObject();
    ltr->dsb_ = dsb_;
    return ltr;
}

void 
dorsalfin::Repeat::copy( const ChromosomeElement &ce ) {
    const Repeat *ltr = dynamic_cast< const Repeat * >( &ce );
    dsb_ = ltr->dsb_;
}

bool
dorsalfin::Repeat::toPool() {
    return ObjectCache< Repeat >::instance()->returnObject( this );
}

void
dorsalfin::Repeat::writeXml( std::ostream &os ) const { 
    os << "<repeat/>\n";
}

