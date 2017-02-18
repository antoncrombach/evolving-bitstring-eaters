//
// Implementation of a grid/population.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "logger.hh"
#include "duo_agent.hh"
#include "agent.hh"
#include "genome.hh"
#include "environment.hh"
#include "population.hh"


bool dorsalfin::Population::shuffle_ = false;
double dorsalfin::Population::threshold_ = 0.0;
int dorsalfin::Population::nr_agent_types_ = 0;
std::string dorsalfin::Population::placement_ = "random";

dorsalfin::Population::Population() 
    : plane_one_(), plane_two_(), 
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(), agent_view_size_( 3 ) {
    scaling_ = new NoScaling();
    selection_ = new RandSelection();
    async_agent_obs_ = 0;
    async_env_change_ = 0;
    async_dsbs_ = 0;
}

dorsalfin::Population::Population( int x, int y, std::vector< Agent* > &vag, 
        ScalingScheme *sca, SelectionScheme *sel ) 
    : plane_one_( boost::extents[ x ][ y ] ), 
      plane_two_( boost::extents[ x ][ y ] ),
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(),  
      scaling_( sca ), selection_( sel ), agent_view_size_( 3 ) {
    // pre: x * y > |vag|
    asLocations( shuffle_locs_ );
    zero( reading );
    zero( writing );
    insert( vag );
    async_agent_obs_ = 0;
    async_env_change_ = 0;
    async_dsbs_ = 0;
}

dorsalfin::Population::Population( int x, int y, std::vector< Agent* > &vag, 
    const std::vector< location > &loc, 
    ScalingScheme *sca, SelectionScheme *sel ) 
    : plane_one_( boost::extents[ x ][ y ] ),
      plane_two_( boost::extents[ x ][ y ] ),
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(),
      scaling_( sca ), selection_( sel ), agent_view_size_( 3 ) {
    // pre: x * y > |vag|
    asLocations( shuffle_locs_ );
    zero( reading );
    zero( writing );
    insertAt( vag, loc );
    async_agent_obs_ = 0;
    async_env_change_ = 0;
    async_dsbs_ = 0;
}

dorsalfin::Population::Population( const Population &pop ) 
    : plane_one_( 
        boost::extents[ 
        pop.write_grid_->shape()[ 0 ] ][ pop.write_grid_->shape()[ 1 ] ] ),
      plane_two_( 
        boost::extents[ 
        pop.write_grid_->shape()[ 0 ] ][ pop.write_grid_->shape()[ 1 ] ] ),
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(), scaling_( 0 ), 
      selection_( 0 ), agent_view_size_( 3 ) {
    // copy deep, all the way down!!
    zero( reading );
    zero( writing );
    // copy agents in map and grid
    for( const_map_ag_iter i = pop.write_agents_.begin(); 
        i != pop.write_agents_.end(); ++i ) {
        // some cloning
        Agent *aux = i->first->clone();
        insertAt( aux, i->second );
    }
    // copy vector of locations
    for( const_loc_iter i = pop.shuffle_locs_.begin(); 
        i != pop.shuffle_locs_.end(); ++i ) {
        shuffle_locs_.push_back( location( *i ) );
    }
    // last but not least...
    scaling_ = pop.scaling_->clone();
    selection_ = pop.selection_->clone();
    model_ = pop.model_;
    async_dsbs_ = 0;
    async_agent_obs_ = 0;
    async_env_change_ = 0;
}

dorsalfin::Population::~Population() {
    // how to deal with the shadow_agents_??
    agents_map aux, bux;
    std::set_intersection( write_agents_.begin(), write_agents_.end(), 
        read_agents_.begin(), read_agents_.end(), 
        std::inserter( aux, aux.begin() ) );
    std::set_symmetric_difference( write_agents_.begin(), write_agents_.end(), 
        read_agents_.begin(), read_agents_.end(), 
        std::inserter( bux, bux.begin() ) );
    // first explicitly delete the agents
    for( const_map_ag_iter i = aux.begin(); i != aux.end(); ++i ) {
        delete i->first;
    }
    for( const_map_ag_iter i = bux.begin(); i != bux.end(); ++i ) {
        delete i->first;
    }
    // and the rest    
    delete scaling_;
    delete selection_;
    if( async_dsbs_ != 0 ) {
        delete async_dsbs_;
    }
    if( async_agent_obs_ != 0 ) {
        delete async_agent_obs_;
    }
    if( async_env_change_ != 0 ) {
        delete async_env_change_;
    }
}

std::vector< dorsalfin::Agent* >
dorsalfin::Population::moore( const location &loc ) {
    // torus border
    int n = read_grid_->shape()[ 0 ];
    int m = read_grid_->shape()[ 1 ];
    // do somethin smart...
    std::vector< Agent* > result;
    int i = (loc.first - 1 + n) % n;
    while( i != (loc.first + 2) % n ) {
        int j = (loc.second - 1 + m) % m;
        while( j != (loc.second + 2) % m ) {
            result.push_back( ( *read_grid_ )[ i ][ j ] );
            j = (j + 1) % m;
        }
        i = (i + 1) % n;
    }
    return result;
}

std::vector< dorsalfin::Agent* >
dorsalfin::Population::moore( const Agent &ag ) {
    return moore( read_agents_[ const_cast< Agent* >( &ag ) ] );
}

dorsalfin::location
dorsalfin::Population::whereIs( const Agent &ag ) {
    return read_agents_[ const_cast< Agent* >( &ag ) ];
}

void
dorsalfin::Population::initialise() {
    // pre: model_ is set
#ifdef DEBUG
    cout << "! agent types: " << nr_agent_types_
         << "\n! agents     : " << write_agents_.size()
         << "\n! placement  : " << placement_ << endl;
#endif
    // do a first shuffle and update of the planes
    swap();
    if( shuffle_ ) shuffle();
    for( map_ag_iter i = write_agents_.begin(); 
        i != write_agents_.end(); ++i ) {
        AgentTag aux = AgentTag( model_->now(), i->second, 0 );
        i->first->myTag( aux );
    }
}

/* hack! */
void
dorsalfin::Population::initAncestorTracing() {
    if( async_agent_obs_ != 0 ) {
        // and init 
        for( map_ag_iter i = write_agents_.begin(); 
            i != write_agents_.end(); ++i ) {
            async_agent_obs_->initialize( i->first );
        }
    }
}
/* end hack */

void 
dorsalfin::Population::evaluate( const Environment &env ) {
    // first dump population in file
    /* hack! */
    if( async_env_change_ != 0 ) {
        async_env_change_->update( this );
    }
    /* end hack */
    for( map_ag_iter i = write_agents_.begin(); 
        i != write_agents_.end(); ++i ) {
        i->first->evaluate( env, i->second );
    }
}

void 
dorsalfin::Population::step() {
    // deterministic synchronous stepping
    // visit every site
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            if( ( *read_grid_ )[ i ][ j ] != 0 ) {
                // occupied spot
                ( *read_grid_ )[ i ][ j ]->step( *this );
                // marked for dying
                if( ( *read_grid_ )[ i ][ j ]->dying() ) {
                    shallowErase( ( *write_grid_ )[ i ][ j ] );
                }
            } else {
                // empty spot
                // get nbh
                location focus( i, j );
                std::vector< Agent* > aux( moore( focus ) );
                // remove non-agents
                aux.erase( std::remove_if( aux.begin(), aux.end(), 
                        IsNotAvailable() ), aux.end() );
                if( aux.size() > 0 ) {
                    // evaluate agents and get their scores
                    std::vector< double > cux;
                    for( ag_iter k = aux.begin(); k < aux.end(); ++k ) {
#ifdef MOORE
                        ( **k ).evaluate( model_->environment(),
                                       whereIs( **k ) );
#else
                        ( **k ).evaluate( model_->environment(), focus );
#endif
                        cux.push_back( ( **k ).score() );
                    }
                    // scale the scores
                    scaling_->scale( cux );
                    // check for sum of fitness
                    double bux = threshold_ -
                        std::accumulate( cux.begin(), cux.end(), 0.0 );
                    if( bux > 0.0 ) {
                        aux.push_back( 0 );
                        cux.push_back( bux );
                    }
                    // select an agent (or a null pointer)
                    Agent *mother = selection_->select( aux, cux );
                    if( mother != 0 ) {
/*                      
                        // log this agent be4 mutations
                        if( async_agent_obs_ != 0 ) {
                            async_agent_obs_->update( mother );
                        }
*/
                        // spawn a sibling
                        Agent *son = mother->sibling();
                        // get the mother tag and set ancestor tags of children
                        AgentTag gux = mother->myTag();
                        mother->parentTag( gux );
                        son->parentTag( gux );
                        // update children's tag
                        if( gux.time == model_->now() ) {
                            // increase index, agent is part of a 'cascade'
                            ++gux.i;
                        }
                        gux.time = model_->now();
                        mother->myTag( gux );
                        son->myTag( AgentTag( model_->now(), focus, 0 ) );
                        // and insert it in the grid
                        insertAt( son, focus );
                        // log after mutations what happened (dsbs)
                        if( async_dsbs_ != 0 ) {
                            DuoAgent hux( mother, son );
                            async_dsbs_->update( &hux );
                        }
                        if( async_agent_obs_ != 0 ) {
                            async_agent_obs_->update( mother );
                            async_agent_obs_->update( son );
                        }
                        // write the produced bitstring into environment
                        model_->environment().update( mother, focus );
                    }
                }
            }
        }
    }
    
    // keep everything consistent
    swap();
    if( shuffle_ ) shuffle();    
}

/*void
dorsalfin::Population::finish() {
    // log all agents to the ancestor tracing observer
    if( async_agent_obs_ != 0 ) {
        for( map_ag_iter i = write_agents_.begin(); 
            i != write_agents_.end(); ++i ) {
            async_agent_obs_->update( i->first );
        }
    }
}*/

void 
dorsalfin::Population::insert( std::vector< Agent* > &ag ) {
    // insert at locations according to placement
    if( placement_ == "patch" ) {
        // assuming 2 different agent types, just sorting them
        std::sort( ag.begin(), ag.end(), SmallerTypeThan() );
        // and need to know how much of each
        ag_iter tt = std::lower_bound( ag.begin(), ag.end(), 
            ag.back(), SmallerTypeThan() );
        std::vector< int > aux;
        aux.push_back( std::distance( ag.begin(), tt ) );
        aux.push_back( std::distance( tt, ag.end() ) );
        // get the patches
        std::vector< std::vector< location > > bux = patchEmptyLocations( aux );
        // only supporting two patches
        insertAt( ag.begin(), tt, bux[ 0 ] );
        insertAt( tt, ag.end(), bux[ 1 ] );
    } else {
        // spread over entire field
        std::vector< location > aux = randEmptyLocations( ag.size() );
        insertAt( ag, aux );
    }
}

void
dorsalfin::Population::insertAt( ag_iter first, ag_iter last, 
        const std::vector< location > &loc ) {
    ag_iter i = first;
    const_loc_iter j = loc.begin();
    while( i != last ) {
        insertAt( *(i++), *(j++) );
    }
}    

void 
dorsalfin::Population::insertAt( 
        std::vector< Agent* > &ag, const std::vector< location > &loc ) {
    // pre: ag.size() == loc.size()
    // pos: memory management of agents is responsibility of this class
    ag_iter i = ag.begin();
    const_loc_iter j = loc.begin();
    while( i != ag.end() ) {
        insertAt( *(i++), *(j++) );
    }
}

void
dorsalfin::Population::insertAt( Agent* ag, const location &loc ) {
    write_agents_[ ag ] = loc;
    ( *write_grid_ )[ loc.first ][ loc.second ] = ag;
}

void 
dorsalfin::Population::erase( std::vector< Agent* > &ag ) {
    std::vector< location > result;
    for( ag_iter i = ag.begin(); i != ag.end(); ++i ) {
        result.push_back( write_agents_[ *i ] );
    }
    eraseAt( result );
}

void
dorsalfin::Population::erase( Agent* ag ) {
    eraseAt( write_agents_[ ag ] );
}

void 
dorsalfin::Population::eraseAt( const std::vector< location > &loc ) {
    for( const_loc_iter i = loc.begin(); i != loc.end(); ++i ) {
        eraseAt( *i );
    }
}

void
dorsalfin::Population::eraseAt( location loc ) {
    write_agents_.erase( ( *write_grid_ )[ loc.first ][ loc.second ] );
    delete ( *write_grid_ )[ loc.first ][ loc.second ];
    ( *write_grid_ )[ loc.first ][ loc.second ] = 0;
}

void
dorsalfin::Population::shallowErase( Agent* ag ) {
    shallowEraseAt( write_agents_[ ag ] );
}

void
dorsalfin::Population::shallowEraseAt( location loc ) {
    write_agents_.erase( ( *write_grid_ )[ loc.first ][ loc.second ] );
    ( *write_grid_ )[ loc.first ][ loc.second ] = 0;
}

std::vector< dorsalfin::location >
dorsalfin::Population::randEmptyLocations( int n ) {
    // find at most n random empty locations
    // note: if necessary, it can be optimised. Yet we're not planning on
    // putting our energy to the optimisation of this method. It's only
    // used once during initialisation of the program FIX
    std::vector< location > result;
    asLocations( result );
    std::random_shuffle( result.begin(), result.end(), rand_range< int > );
    loc_iter i( result.begin() );
    loc_iter j( boost::prior( result.end() ) );
    while( i != j ) {
        if( ( *write_grid_ )[ i->first ][ i->second ] == 0 ) {
            ++i;
        } else {
            std::iter_swap( i, j );
            --j;
        }
    }

    int k = std::distance( result.begin(), i );
    if( k < n ) n = k;
    result.erase( result.begin() + n, result.end() );
    return result;
}

std::vector< std::vector< dorsalfin::location > >
dorsalfin::Population::patchEmptyLocations( std::vector< int > n ) {
    // find n random locations near each other
    // no general algorithm... only nr_agent_types = {1, 2} supported
    std::vector< std::vector< location > > result;
    std::vector< location > aux;
    asLocations( aux );
    std::random_shuffle( aux.begin(), aux.end(), rand_range< int > );
    int bux = static_cast< int >( write_grid_->shape()[ 0 ] ) / 2;
    // do it
    for( std::vector< int >::iterator ii = n.begin(); ii != n.end(); ++ii ) {
        loc_iter jj = aux.begin();
        result.push_back( std::vector< location >() );
        int k = 0; 
        while( k != *ii ) {
            if( jj->first >= std::distance( n.begin(), ii ) * bux and
                jj->first < ( std::distance( n.begin(), ii ) + 1 ) * bux ) {
                result.back().push_back( *jj );
                ++k;
            }
            // no check if jj goes past end!!
            ++jj;
        }
    }
    return result;
}

void
dorsalfin::Population::asLocations( std::vector< location > &ll ) {
    // grids are of same size
    ll.clear();
    ll.reserve( read_grid_->shape()[ 0 ] * read_grid_->shape()[ 1 ] );
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            ll.push_back( std::make_pair( i, j ) );
        }
    }
}

void
dorsalfin::Population::zero( grid_type tt ) {
    agents_grid *g = write_grid_;
    if( tt == reading ) {
        g = read_grid_;
    }
    
    for( uint i = 0; i < g->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < g->shape()[ 1 ]; ++j ) {
            ( *g )[ i ][ j ] = static_cast< Agent* >( 0 );
        }
    }
}

void
dorsalfin::Population::swap() {
    // swap grids and mapping, shadow is readable now...
    std::swap( write_grid_, read_grid_ );
    write_agents_.swap( read_agents_ );
    
    // keep data consistent
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            if( (* read_grid_ )[ i ][ j ] != 0 and 
                ( *write_grid_ )[ i ][ j ] == 0 ) {
                // insertion happened last time
                insertAt( ( *read_grid_ )[ i ][ j ], std::make_pair( i, j ) ); 
            } else if( (* read_grid_ )[ i ][ j ] == 0 and
                ( *write_grid_ )[ i ][ j ] != 0 ) {
                // deletion happened last time
                /* ugly hack, but correct */
                if( async_agent_obs_ != 0 ) {
                    async_agent_obs_->update( ( *write_grid_ )[ i ][ j ] );
                }
                /* end hack */
                eraseAt( std::make_pair( i, j ) );
            }
        }
    }
}

void
dorsalfin::Population::swap( grid_type tt, const location &a, const location &b ) {
    agents_grid *g = write_grid_;
    agents_map *m = &write_agents_;
    if( tt == reading ) {
        g = read_grid_;
        m = &read_agents_;
    }
    // first swap in map
    if( ( *g )[ a.first ][ a.second ] != 0 ) {
        ( *m )[ ( *g )[ a.first ][ a.second ] ] = b;
    }
    if( ( *g )[ b.first ][ b.second ] != 0 ) {
        ( *m )[ ( *g )[ b.first ][ b.second ] ] = a;
    }
    // now swap in grid
    std::swap( ( *g )[ a.first ][ a.second ], ( *g )[ b.first ][ b.second ] );
}

void
dorsalfin::Population::shuffle() {
    // assuming grids are consistent
    std::random_shuffle( shuffle_locs_.begin(), shuffle_locs_.end(),  
        rand_range< int > );
    const_loc_iter aux = shuffle_locs_.begin();
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            swap( reading, *aux, std::make_pair( i, j ) );
            swap( writing, *aux, std::make_pair( i, j ) );
            ++aux;
        }
    }
}

void
dorsalfin::Population::attach1( AsyncLogObserver *l ) {
    async_agent_obs_ = l;
}

void
dorsalfin::Population::attach2( AsyncLogObserver *l ) {
    async_env_change_ = l;
}

void
dorsalfin::Population::attach3( AsyncLogObserver *l ) {
    async_dsbs_ = l;
}

void
dorsalfin::Population::closeAll() {
    if( async_agent_obs_ != 0 ) {
        async_agent_obs_->closeLog();
    }
    if( async_env_change_ != 0 ) {
        async_env_change_->closeLog();
    }
    if( async_dsbs_ != 0 ) {
        async_dsbs_->closeLog();
    }
}

void
dorsalfin::Population::threshold( double f ) 
{ threshold_ = f; }
    
double
dorsalfin::Population::threshold()
{ return threshold_; }

void
dorsalfin::Population::nrAgentTypes( int t )
{ nr_agent_types_ = t; }

int
dorsalfin::Population::nrAgentTypes()
{ return nr_agent_types_; }

void
dorsalfin::Population::placement( const std::string &s )
{ placement_ = s; }

std::string
dorsalfin::Population::placement()
{ return placement_; }

void
dorsalfin::Population::shuffling( bool t )
{ shuffle_ = t; }

bool
dorsalfin::Population::shuffling()
{ return shuffle_; }

bool
dorsalfin::Population::hasEveryAgentType() const {
    std::vector< uint > aux( nr_agent_types_, 0 );
    const_map_ag_iter i = write_agents_.begin();
    const_map_ag_iter j = write_agents_.end();
    while( i != j ) {
        if( std::accumulate( aux.begin(), aux.end(), 0 ) != nr_agent_types_ ) {
            aux[ i->first->type() ] = 1;
            ++i;
        } else {
            j = i;
        }
    }
    return std::accumulate( aux.begin(), aux.end(), 0 ) == nr_agent_types_;
}

void
dorsalfin::Population::writeXml( std::ostream &os ) const {
    for( const_map_ag_iter i = write_agents_.begin(); 
        i != write_agents_.end(); ++i ) {
        // output the genome agents
        os << *( i->first );
    }
}
