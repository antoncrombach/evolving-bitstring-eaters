//
// Implementation of functions of any network agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "genreg_agent.hh"
#include "environment.hh"
#include "population.hh"


double dorsalfin::GenRegAgent::death_rate_ = 0.0;

int dorsalfin::GenRegAgent::max_dist_ = 0;
boost::dynamic_bitset<> dorsalfin::GenRegAgent::init_state_( 0 );

int dorsalfin::GenRegAgent::penalty_genome_ = 0;
int dorsalfin::GenRegAgent::penalty_tposons_ = 0;
double dorsalfin::GenRegAgent::penalty_genome_rate_ = 0;
double dorsalfin::GenRegAgent::penalty_tposons_rate_ = 0;


dorsalfin::GenRegAgent::GenRegAgent() 
    : Agent( 0 ), distance_( 0 ), 
      genome_( 0 ), network_( 0 ), read_from_file_( false ) {}

dorsalfin::GenRegAgent::GenRegAgent( uint tag ) 
    : Agent( tag ), distance_( 0 ), 
      genome_( 0 ), network_( 0 ), read_from_file_( false ) {}
      
dorsalfin::GenRegAgent::GenRegAgent( uint tag, Genome *g ) 
    : Agent( tag ), distance_( 0 ), 
      genome_( g ), network_( new Network() ), 
      read_from_file_( false ) {}
    
dorsalfin::GenRegAgent::GenRegAgent( const GenRegAgent &ag ) 
    : Agent( 0 ), genome_( 0 ), network_( 0 ) {
    GenRegAgent::copy( ag );
}

dorsalfin::GenRegAgent::~GenRegAgent() {
    if( genome_ != 0 ) delete genome_;
    if( network_ != 0 ) delete network_;
}

void 
dorsalfin::GenRegAgent::copy( const Agent &ag ) {
    const GenRegAgent &sag = dynamic_cast< const GenRegAgent & >( ag );
    Agent::copy( ag );
    distance_ = sag.distance_;
    parent_.size_ = sag.parent_.size_;
    parent_.distance_ = sag.parent_.distance_;
    genome_ = sag.genome_->clone();
    // note: graph needs to be rebuilt!!
    network_ = sag.network_->clone();
    read_from_file_ = sag.read_from_file_;
}

void
dorsalfin::GenRegAgent::initialise() {
    if( not read_from_file_ ) {
        genome_->setState( init_state_ );
        genome_->initialise();
    }
    // Make sure caching is ok
    // genome stuff
    genome_->nrRetroposons();
    genome_->length();
    // network stuff.. not necessary 
}

double
dorsalfin::GenRegAgent::score() const {
    double result = static_cast< double >( distance_ );
    // now see if we have any penalties to add
    int aux = genome_->length() - penalty_genome_;
    int bux = genome_->nrRetroposons() - penalty_tposons_;
    if( aux > 0 ) {
        result += aux * penalty_genome_rate_;
    }
    if( bux > 0 ) {
        result += bux * penalty_tposons_rate_;
    }
    // we've got the raw score now and calculate the fitness score
    if( result < max_dist_ ) {
        return 1.0 - ( result / static_cast< double >( max_dist_ ) );
    } else {
        return 0.0;
    }
}

void
dorsalfin::GenRegAgent::deathRate( double f )
{ death_rate_ = f; }

double
dorsalfin::GenRegAgent::deathRate()
{ return death_rate_; }

void
dorsalfin::GenRegAgent::maxDistance( int f ) 
{ max_dist_ = f; }

int
dorsalfin::GenRegAgent::maxDistance() 
{ return max_dist_; }

void
dorsalfin::GenRegAgent::maxGenomeSize( int f )
{ penalty_genome_ = f; }

void
dorsalfin::GenRegAgent::maxTposons( int f )
{ penalty_tposons_ = f; }

void
dorsalfin::GenRegAgent::genomePenaltyRate( double f )
{ penalty_genome_rate_ = f; }

void
dorsalfin::GenRegAgent::retroposonPenaltyRate( double f )
{ penalty_tposons_rate_ = f; }

void 
dorsalfin::GenRegAgent::initialState( const boost::dynamic_bitset<> &f )
{ init_state_ = f; }

