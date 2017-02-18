//
// Implementation of a special chromosome element, the interaction.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "interaction.hh"
#include "pool.hh"

template<> dorsalfin::ObjectCache< dorsalfin::Interaction >* 
dorsalfin::ObjectCache< dorsalfin::Interaction >::instance_ = 0;

dorsalfin::Interaction::Interaction() 
    : ChromosomeElement(), weight_( -1 ), tag_( 0 ) {}

dorsalfin::Interaction::Interaction( int t, int w )
    : ChromosomeElement(), weight_( w ), tag_( t ) {}
    
dorsalfin::Interaction::Interaction( const Interaction &tp ) 
    : ChromosomeElement( tp ) {
    Interaction::copy( tp );
}

dorsalfin::ChromosomeElement* 
dorsalfin::Interaction::clone() const {
    //return new Interaction( *this );
    Interaction *ia = ObjectCache< Interaction >::instance()->borrowObject();
    ia->weight_ = weight_;
    ia->tag_ = tag_;
    return ia;
}

void 
dorsalfin::Interaction::copy( const ChromosomeElement &ce ) {
    const Interaction *ia = dynamic_cast< const Interaction * >( &ce );
    weight_ = ia->weight_;
    tag_ = ia->tag_;
}

bool
dorsalfin::Interaction::toPool() {
    return ObjectCache< Interaction >::instance()->returnObject( this );
}

void
dorsalfin::Interaction::writeXml( std::ostream &os ) const { 
    os << "<ia t=\"" << tag_ << "\" w=\"" << weight_ << "\"/>\n";
}

