//
// Implementation of a simple agent (not complete).
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "simple_agent.hh"
#include "population.hh"


double dorsalfin::SimpleAgent::birth_rate_ = 0.0;
double dorsalfin::SimpleAgent::death_rate_ = 0.0;


dorsalfin::SimpleAgent::SimpleAgent() : Agent( 0 ), score_( birth_rate_ ) {}

dorsalfin::SimpleAgent::SimpleAgent( const SimpleAgent &ag ) : Agent( 0 ) {
    SimpleAgent::copy( ag );
}

dorsalfin::Agent* 
dorsalfin::SimpleAgent::clone() const {
    return new SimpleAgent( *this );
}

void 
dorsalfin::SimpleAgent::copy( const Agent &ag ) {
    const SimpleAgent &sag = dynamic_cast< const SimpleAgent & >( ag );
    Agent::copy( ag );
    score_ = sag.score_;
}

void 
dorsalfin::SimpleAgent::step( Population &pop ) {
    double rr = uniform();
    if( rr < death_rate_ ) {
        dying_ = true;
    }
}

dorsalfin::Agent*
dorsalfin::SimpleAgent::sibling() {
    Agent* result;
    result = this->clone();
    //SimpleAgent *aux = dynamic_cast< SimpleAgent* >( result );
    return result;
}

void 
dorsalfin::SimpleAgent::writeXml( std::ostream &os ) const {
    os << "<agent birth=\"" << me_.time << "\" x=\"" << me_.loc.first;
    os << "\" y=\"" << me_.loc.second << "\" i=\"" << me_.i << "\"";
    os << ">\n<class>SimpleAgent</class>\n</agent>\n";
}

void
dorsalfin::SimpleAgent::birthRate( double f )
{ birth_rate_ = f; }

double
dorsalfin::SimpleAgent::birthRate()
{ return birth_rate_; }

void
dorsalfin::SimpleAgent::deathRate( double f )
{ death_rate_ = f; }

double
dorsalfin::SimpleAgent::deathRate()
{ return death_rate_; }

