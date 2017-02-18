//
// Implementation of a special gene, the Retroposon.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "retroposon.hh"

template<> dorsalfin::ObjectCache< dorsalfin::Retroposon >* 
dorsalfin::ObjectCache< dorsalfin::Retroposon >::instance_ = 0;


dorsalfin::ChromosomeElement* 
dorsalfin::Retroposon::clone() const {
    //return new Retroposon( *this );
    Retroposon *rp = ObjectCache< Retroposon >::instance()->borrowObject();
    rp->copy( *this );
    return rp;
}

void 
dorsalfin::Retroposon::copy( const ChromosomeElement &ce ) {
    const Retroposon *tp = dynamic_cast< const Retroposon * >( &ce );
    tag_ = tp->tag_;
}

bool
dorsalfin::Retroposon::toPool() {
    return ObjectCache< Retroposon >::instance()->returnObject( this );
}

void
dorsalfin::Retroposon::writeXml( std::ostream &os ) const {
    os << "<tposon id=\"" + boost::lexical_cast< std::string >( tag_ ) + "\"/>";
}
