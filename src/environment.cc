//
// Implementation of bitset environment
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "environment.hh"
#include "distribution.hh"

dorsalfin::Environment::Environment( uint s ) 
    : generator_env_( s ), uniform_env_( generator_env_ ), model_( 0 ) {}

void
dorsalfin::Environment::finish() {
    notify();
}

//
// Bitset environment
//
double dorsalfin::BitEnvironment::diffusion_ = 0.0;

dorsalfin::BitEnvironment::BitEnvironment( int x, int y, uint s ) 
    : Environment( s ), rw_grid_( boost::extents[ x ][ y ] ), 
      size_( 0 ), period_( 0 ), bits_() {}
   
void
dorsalfin::BitEnvironment::initialise( 
    const std::vector< boost::dynamic_bitset<> > &b ) {
    std::copy( b.begin(), b.end(), std::back_inserter( bits_ ) );   
    size_ = b.front().size();
}

void
dorsalfin::BitEnvironment::moore( 
    std::vector< boost::dynamic_bitset<> > &result, const location &loc) const {
    // torus border
    int n = rw_grid_.shape()[ 0 ];
    int m = rw_grid_.shape()[ 1 ];
    // do somethin smart...
    int i = (loc.first - 1 + n) % n;
    while( i != (loc.first + 2) % n ) {
        int j = (loc.second - 1 + m) % m;
        while( j != (loc.second + 2) % m ) {
            result.push_back( rw_grid_[ i ][ j ] );
            j = (j + 1) % m;
        }
        i = (i + 1) % n;
    }
}

void
dorsalfin::BitEnvironment::diffusion() {
    // issue warning
    uint x = rw_grid_.shape()[ 0 ];
    uint y = rw_grid_.shape()[ 1 ];
    if( x & 1 or y & 1 ) {
        std::cerr << "Warning: cannot do diffusion, uneven x,y" << std::endl;
        return;
    }
    // Margolus diffusion
    for( uint i = period_; i < x; i += 2 ) {
        for( uint j = period_; j < y; j += 2 ) {
            rotate( i, j, x, y );
        }
    }
    period_ = 1 - period_;
}

inline void
dorsalfin::BitEnvironment::rotate( uint i, uint j, uint x, uint y ) {
    // pre: left upper corner is given by (i,j)
    uint ib = ( i + 1 + x ) % x;
    uint jc = ( j + 1 + y ) % y;
    
    if( uniform_env_() < 0.5 ) {
        rw_grid_[ i ][ j ].swap( rw_grid_[ ib ][ j ] );
        rw_grid_[ ib ][ j ].swap( rw_grid_[ ib ][ jc ] );
        rw_grid_[ ib ][ jc ].swap( rw_grid_[ i ][ jc ] );
    } else {
        rw_grid_[ i ][ j ].swap( rw_grid_[ i ][ jc ] );
        rw_grid_[ i ][ jc ].swap( rw_grid_[ ib ][ jc ] );
        rw_grid_[ ib ][ jc ].swap( rw_grid_[ ib ][ j ] );
    }
}

void
dorsalfin::BitEnvironment::setDiffusion( double l ) {
    diffusion_ = l;
}

void
dorsalfin::BitEnvironment::readGrid( std::ifstream &infile ) {
    for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
            std::string aux;
            infile >> aux;
            rw_grid_[ i ][ j ] = boost::dynamic_bitset<>( aux.erase( 0, 2 ) );
        }
    }
}

//
// Noisy bitset environment
//
dorsalfin::NoisyBitEnvironment::NoisyBitEnvironment()
    : BitEnvironment() {
    generateLocations();
}

dorsalfin::NoisyBitEnvironment::NoisyBitEnvironment( int x, int y, 
    uint s, double r ) 
    : BitEnvironment( x, y, s ), rate_( r ) {
    generateLocations();
}

void
dorsalfin::NoisyBitEnvironment::initialise( 
    const std::vector< boost::dynamic_bitset<> > &b ) {
    BitEnvironment::initialise( b );
    if( b.size() == 1 ) {
        // only 1, homogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = bits_.front();
            }
        }
    } else {
        // multiple, well mixed heterogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = 
                    *( random_element( 
                        bits_.begin(), bits_.end(), rand_range< int > ));
            }
        }
    }
}

void
dorsalfin::NoisyBitEnvironment::fluctuate( long time ) {
    // how many sites to update?
    if( not close_to( rate_, 0.0 ) ) {
        uint n = Binomial< uniform_gen_type >()( cache_locs_.size(), rate_,
            uniform_env_ );
        // get some sites
        std::vector< location > sample( n );
        __gnu_cxx::random_sample( cache_locs_.begin(), cache_locs_.end(),
            sample.begin(), sample.end(), rand_range< int > );
        for( std::vector< location >::iterator i = sample.begin(); 
            i != sample.end(); ++i ) {
            dorsalfin::rotate( rw_grid_[ i->first ][ i->second ], 1 );
        }
    }
    // and diffuse
    if( diffusion_ > 1.0 ) {
        for( int i = 0; i != static_cast< int >( diffusion_ ); ++i ) {
            // add i to make the even/uneven thing work 
            diffusion();
        }
    } else {
        if( not close_to( diffusion_, 0.0 ) ) {
            if( time % static_cast< int >( 1 / diffusion_ ) == 0 ) {
                diffusion();
            }
        }
    }
}

void
dorsalfin::NoisyBitEnvironment::update( const Agent *ag, 
    const location &loc ) {
    const EnvAgent *eag = dynamic_cast< const EnvAgent* >( ag );
    // check if we return quickly as agent is not fit enough
    if( eag->nrPredictedBits() < 1 ) return;
    // get input pattern from location and agent
    boost::dynamic_bitset<> in;
    getBits( in, loc );
    // rotate what we understood
    dorsalfin::rotate( in, eag->nrPredictedBits() );
    setBits( in, loc );
}

void
dorsalfin::NoisyBitEnvironment::generateLocations() {
    for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
            cache_locs_.push_back( std::make_pair( i, j ) );
        }
    }
}

//
// NoRotate bitset environment
//
dorsalfin::NoRotateBitEnvironment::NoRotateBitEnvironment()
    : BitEnvironment() {}

dorsalfin::NoRotateBitEnvironment::NoRotateBitEnvironment( int x, int y, 
    uint s, double r ) : BitEnvironment( x, y, s ), influx_( r ) {}

void
dorsalfin::NoRotateBitEnvironment::initialise( 
    const std::vector< boost::dynamic_bitset<> > &b ) {
    BitEnvironment::initialise( b );
    if( b.size() == 1 ) {
        // only 1, homogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = b.front();
            }
        }
    } else {
        // multiple, well mixed heterogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = 
                    *( random_element( b.begin(), b.end(), rand_range< int > ));
            }
        }
    }
}
    
void
dorsalfin::NoRotateBitEnvironment::fluctuate( long time ) {
    // influx of new bitstrings
    if( not close_to( influx_, 0.0 ) ) {
        // probalistic adding of bitsets
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                if( rw_grid_[ i ][ j ].empty() and uniform_env_() < influx_ ) {
                     rw_grid_[ i ][ j ] = *( random_element( bits_.begin(),
                         bits_.end(), rand_range< int > ) );
                    
                }
            }
        }
    } else {
        // replenish everything
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                if( rw_grid_[ i ][ j ].empty() ) {
                     rw_grid_[ i ][ j ] = *( random_element( bits_.begin(),
                         bits_.end(), rand_range< int > ) );
                }
            }
        }
    }
    // and diffuse
    if( diffusion_ > 1.0 ) {
        for( int i = 0; i != static_cast< int >( diffusion_ ); ++i ) {
            // add i to make the even/uneven thing work 
            diffusion();
        }
    } else {
        if( not close_to( diffusion_, 0.0 ) ) {
            if( time % static_cast< int >( 1 / diffusion_ ) == 0 ) {
                diffusion();
            }
        }
    }
}

void
dorsalfin::NoRotateBitEnvironment::update( const Agent *ag, 
    const location &loc ) {
    const EnvAgent *eag = dynamic_cast< const EnvAgent* >( ag );
    // check if we return quickly as agent is not fit enough
    if( eag->nrPredictedBits() < 1 ) return;
    // get input pattern from location and agent
    boost::dynamic_bitset<> in;
    getBits( in, loc );
    // remove what we understood
    in >>= eag->nrPredictedBits();
    in.resize( in.size() - eag->nrPredictedBits() );
    setBits( in, loc );
}


//
// Subset bitset environment
//
dorsalfin::SubsetBitEnvironment::SubsetBitEnvironment()
    : BitEnvironment() {}

dorsalfin::SubsetBitEnvironment::SubsetBitEnvironment( int x, int y, uint s ) 
    : BitEnvironment( x, y, s ) {}

void
dorsalfin::SubsetBitEnvironment::initialise( 
    const std::vector< boost::dynamic_bitset<> > &b ) {
    BitEnvironment::initialise( b );
    if( b.size() == 1 ) {
        // only 1, homogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = b.front();
            }
        }
    } else {
        // multiple, well mixed heterogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = 
                    *( random_element( b.begin(), b.end(), rand_range< int > ));
            }
        }
    }
}

void
dorsalfin::SubsetBitEnvironment::fluctuate( long time ) {
    if( diffusion_ > 1.0 ) {
        for( int i = 0; i != static_cast< int >( diffusion_ ); ++i ) {
            // add i to make the even/uneven thing work 
            diffusion();
        }
    } else {
        if( not close_to( diffusion_, 0.0 ) ) {
            if( time % static_cast< int >( 1 / diffusion_ ) == 0 ) {
                diffusion();
            }
        }
    }
}

//
// OneBite bitset environment
//
dorsalfin::OneBiteBitEnvironment::OneBiteBitEnvironment()
    : BitEnvironment() {}

dorsalfin::OneBiteBitEnvironment::OneBiteBitEnvironment( int x, int y, 
    uint s, double r ) : BitEnvironment( x, y, s ), influx_( r ) {}

void
dorsalfin::OneBiteBitEnvironment::initialise( 
    const std::vector< boost::dynamic_bitset<> > &b ) {
    BitEnvironment::initialise( b );
    if( b.size() == 1 ) {
        // only 1, homogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = b.front();
            }
        }
    } else {
        // multiple, well mixed heterogeneous grid
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                rw_grid_[ i ][ j ] = 
                    *( random_element( b.begin(), b.end(), rand_range< int > ));
            }
        }
    }
}

void
dorsalfin::OneBiteBitEnvironment::fluctuate( long time ) {
    // influx of new bitstrings
    if( not close_to( influx_, 0.0 ) ) {
        // probalistic adding of bitsets
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                if( rw_grid_[ i ][ j ].empty() and uniform_env_() < influx_ ) {
                     rw_grid_[ i ][ j ] = *( random_element( bits_.begin(),
                         bits_.end(), rand_range< int > ) );
                    
                }
            }
        }
    } else {
        // replenish everything
        for( uint i = 0; i < rw_grid_.shape()[ 0 ]; ++i ) {
            for( uint j = 0; j < rw_grid_.shape()[ 1 ]; ++j ) {
                if( rw_grid_[ i ][ j ].empty() ) {
                     rw_grid_[ i ][ j ] = *( random_element( bits_.begin(),
                         bits_.end(), rand_range< int > ) );
                }
            }
        }
    }
    if( diffusion_ > 1.0 ) {
        for( int i = 0; i != static_cast< int >( diffusion_ ); ++i ) {
            // add i to make the even/uneven thing work 
            diffusion();
        }
    } else {
        if( not close_to( diffusion_, 0.0 ) ) {
            if( time % static_cast< int >( 1 / diffusion_ ) == 0 ) {
                diffusion();
            }
        }
    }
}

void
dorsalfin::OneBiteBitEnvironment::update( const Agent *ag, 
    const location &loc ) {
    const EnvAgent *eag = dynamic_cast< const EnvAgent* >( ag );
    // check if we return quickly as agent is not fit enough
    if( eag->nrPredictedBits() < 1 ) return;
    // set empty input pattern from location and agent
    boost::dynamic_bitset<> in;
    setBits( in, loc );
}


// below we have some old update mechanisms that did not work properly
/*void
dorsalfin::NoisyBitEnvironment::update( const Agent *ag, 
    const location &loc ) {
    // prob of writing
    const double write_ = 0.01;
    const EnvAgent *eag = dynamic_cast< const EnvAgent* >( ag );
    // check if we return quickly as agent is not fit enough
    if( eag->nrPredictedBits() == 0 ) return;
    // let's go!
    if( uniform_env_() < write_ ) {
        // get input pattern from location and agent
        boost::dynamic_bitset<> in;
        getBits( in, loc );
        boost::dynamic_bitset<> out( eag->output() );
        out.resize( size_ );
        // append and shift part of in
        out <<= ( size_ - eag->nrPredictedBits() );
        in >>= eag->nrPredictedBits();
        in |= out;
        setBits( in, loc );
    }
}*/

/*void
dorsalfin::NoisyBitEnvironment::update( const Agent *ag, 
    const location &loc ) {
    // prob of writing
    const double write_ = 0.01;
    const EnvAgent *eag = dynamic_cast< const EnvAgent* >( ag );
    // check if we return quickly as agent is not fit enough
    if( eag->distance() == GenRegAgent::maxDistance() ) return;
    // let's go!
    if( uniform_env_() < write_ ) {
        // get input pattern from location and agent
        boost::dynamic_bitset<> in;
        getBits( in, loc );
        // resize out
        boost::dynamic_bitset<> out( eag->output() );
        int d = out.size() - in.size();
        if( d > 0 ) {
            out >>= d;
            out.resize( in.size() );
        }
        // overwrite part of in
        boost::dynamic_bitset<> mask( in.size(), 0ul );
        for( uint i = 0; i != ( in.size() - eag->distance() ); ++i ) {
            mask.set( i );
        }
        in = ( in & ~mask ) | ( out & mask );
        setBits( in, loc );
    }
}*/
