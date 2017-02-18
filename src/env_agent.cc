//
// Implementation of a network agent with environmental interactions.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "env_agent.hh"
#include "environment.hh"
#include "population.hh"


dorsalfin::EnvAgent::EnvAgent() : GenRegAgent() {}

dorsalfin::EnvAgent::EnvAgent( uint tag, Genome *g )
    : GenRegAgent( tag, g ), reproduced_( false ) {}

dorsalfin::EnvAgent::EnvAgent( const EnvAgent &ag ) 
    : GenRegAgent( ag ), reproduced_( false ) {
    EnvAgent::copy( ag );
}

dorsalfin::EnvAgent::~EnvAgent() {}

dorsalfin::Agent* 
dorsalfin::EnvAgent::clone() const {
    return new EnvAgent( *this );
}

void 
dorsalfin::EnvAgent::copy( const Agent &ag ) {}

void 
dorsalfin::EnvAgent::step( Population &pop ) {
    // are we dead?
    double rr = uniform();
    if( rr < death_rate_ or distance_ < 0 ) {
        dying_ = true;
    }
}

dorsalfin::Agent*
dorsalfin::EnvAgent::sibling() {
    // change my parent info, coz after mitosis I'm not parent anymore
    parent_.size_ = genome_->length();
    parent_.distance_ = distance_;
    // first mitosis: now genome changes
    genome_->duplicate();
    genome_->mutate();
    Genome *sister_genome = genome_->split();
    // building sister net agent
    EnvAgent *sister = new EnvAgent( type_, sister_genome );
    // update parent info in sister
    sister->parent_.size_ = parent_.size_;
    sister->parent_.distance_ = parent_.distance_;
    // and become reborn
    reproduced_ = true;
    return sister;
}

void
dorsalfin::EnvAgent::evaluate( const Environment &env, const location &loc ) {
    // Not anymore just born, trying to reproduce
    reproduced_ = false;
    just_born_ = false;
    // Not just any environment
    const BitEnvironment &be = 
        dynamic_cast< const BitEnvironment & >( env );
    // and evaluate agent
    if( network_->isEmpty() or genome_->hasMutation() ) {
        network_->build( *genome_ );
    }
    // successful build?
    if( network_->hasOutputNode() ) {
        boost::dynamic_bitset<> target;
        be.getBits( target, loc );
        // set states of genes 
        // Note: now they are first all set to zero!
        network_->resetState();
        network_->setInput( target );
        network_->setInputState( target );
        // do network propagation
        network_->propagate();
        // and evaluate fitness
        distance_ = network_->output().size();
    } else {
        distance_ = -1;
    }
}

double
dorsalfin::EnvAgent::score() const {
    // return quickly if not evaluated yet
    if( just_born_ ) return 0.0;
    // and go
    double result = static_cast< double >( distance_ );
    // now see if we have any penalties to add to the distance
    int aux = genome_->length() - penalty_genome_;
    int bux = genome_->nrRetroposons() - penalty_tposons_;
    if( aux > 0 ) {
        result -= aux * penalty_genome_rate_;
    }
    if( bux > 0 ) {
        result -= bux * penalty_tposons_rate_;
    }
    // and do some exponential stuff, which needs linear scaling...
    return std::max( 0.0, result );
}

void 
dorsalfin::EnvAgent::writeXml( std::ostream &os ) const {
    // pre: genome_ != 0
    os << "<agent birth=\"" << me_.time 
       << "\" x=\"" << me_.loc.first << "\" y=\"" << me_.loc.second
       << "\" i=\"" << me_.i << "\">\n";
    os << "<parent birth=\"" << ancestor_.time
       << "\" x=\"" << ancestor_.loc.first << "\" y=\"" << ancestor_.loc.second
       << "\" i=\"" << ancestor_.i << "\"/>\n<class>EnvAgent</class>\n";
    os << "<score>" << score() << "</score>\n"
       << "<dist>" << distance_ << "</dist>\n" << *genome_; 
    if( reproduced_ ) {
        os << "<out>" << network_->output() << "</out>\n";
    } else {
        os << "<out>none</out>\n";
    }
    os << "</agent>\n";
}


