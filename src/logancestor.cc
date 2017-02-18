//
// Implementation of a few observers that dump their data in files.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "logancestor.hh"

//#include "environment.hh"
#include "population.hh"
#include "genreg_agent.hh"
#include "genome.hh"
//#include "chromosome.hh"
//#include "retroposon.hh"
//#include "gene.hh"

//#include "duo_agent.hh"


//
// More complicated csv observer for ancestor tracing
//
// Note: have to make an educated guess for the max allowed size,
// should depend on grid size...
uint dorsalfin::LogCsvAncestors::max_forest_size_ = 100000;

dorsalfin::LogCsvAncestors::LogCsvAncestors( 
    std::string fname ) : AsyncLogObserver(), nr_roots_( 0 ), forest_(),
    leafs_(), mothers_(), vertex_nr_( 0 ), edge_nr_( 0 ) {
    openLog( fname );
    writeHeader();
    initForest();
}

void
dorsalfin::LogCsvAncestors::initialize( Subject *s ) {
    // init agent by putting its id in the leafs and forest
    Agent *ag = dynamic_cast< Agent * >( s );
    AgentTag ch = ag->myTag();
    AgentTag mo;
    
    leafs_.insert( std::make_pair( ch, boost::add_vertex( forest_ ) ) );
    boost::put( boost::vertex_index, forest_, leafs_[ ch ], vertex_nr_++ );
    boost::put( boost::vertex_bundle, forest_, leafs_[ ch ], ch );
    // \c mo is null object and present in forest after initForest
    anc_edge aux1 = boost::add_edge( leafs_[ ch ], mothers_[ mo ], 
        forest_ ).first;
    boost::put( boost::edge_index, forest_, aux1, edge_nr_++ );
    ++nr_roots_;
}

void 
dorsalfin::LogCsvAncestors::update( Subject *s ) {
    Agent *ag = dynamic_cast< Agent * >( s );
    if( nr_roots_ == 10 ) {
        // renew buffer, write old contents
        writeForest();
        renewForest();
    } 
    if( boost::num_vertices( forest_ ) > max_forest_size_ ) {
        // if forest becomes really large...
        writeForest();
        renewForest();
    }

    // what to do?
    if( leafs_.find( ag->myTag() ) != leafs_.end() ) {
        death( ag );
    } else {
        birth( ag );
    }
}

void
dorsalfin::LogCsvAncestors::birth( Agent *ag ) {
    AgentTag ch = ag->myTag();
    AgentTag mo = ag->parentTag();

    std::map< AgentTag, anc_vertex >::iterator aux = leafs_.find( mo );
    if( aux != leafs_.end() ) {
        // not a mother yet, check grandma too (if there is any...)
        AgentTag g;
        adj_iter vg, end;
        boost::tie( vg, end ) = boost::adjacent_vertices( aux->second,forest_ );
        // assert vg != end, ie there is a granny
        if( vg != end ) {
            g = boost::get( boost::vertex_bundle, forest_, *vg );
            invadj_iter i, j, k;
            boost::tie( i, j ) = boost::inv_adjacent_vertices( *vg, forest_ );
            k = j;
            while( i != j ) {
                std::map< AgentTag, anc_vertex >::iterator bux =
                    leafs_.find( 
                        boost::get( boost::vertex_bundle, forest_, *i ) );
                if( bux == aux or bux == leafs_.end() ) { 
                    ++i;
                } else {
                    j = i;
                }
            }
            if( j == k ) {
                // no more leafs to be found, grandma shouldn't be a mom
                mothers_.erase( g );
            }
        } else {
            std::cout << "No granny" << std::endl;
        }
        // add mo to mother and remove it from leafs
        mothers_.insert( *aux );
        leafs_.erase( aux );
    } else if( mothers_.find( mo ) == mothers_.end() ) {
        // only happens if inbetween two births the forest has been written 
        AgentTag rt;
        mothers_.insert( std::make_pair( mo, boost::add_vertex( forest_ ) ) );
        boost::put( boost::vertex_index, forest_, mothers_[ mo ],   
            vertex_nr_++ );
        boost::put( boost::vertex_bundle, forest_, mothers_[ mo ], mo );
        anc_edge aux1 = boost::add_edge( mothers_[ mo ], mothers_[ rt ], 
            forest_ ).first;
        boost::put( boost::edge_index, forest_, aux1, edge_nr_++ );
    }

    // and add the child (mo already in mothers)
    leafs_.insert( std::make_pair( ch, boost::add_vertex( forest_ ) ) );
    boost::put( boost::vertex_index, forest_, leafs_[ ch ], vertex_nr_++ );
    boost::put( boost::vertex_bundle, forest_, leafs_[ ch ], ch );
    anc_edge aux1 = boost::add_edge( leafs_[ ch ], mothers_[ mo ], 
        forest_ ).first;
    boost::put( boost::edge_index, forest_, aux1, edge_nr_++ );
}

void
dorsalfin::LogCsvAncestors::death( Agent *ag ) {
    AgentTag ch = ag->myTag();
    AgentTag mo = ag->parentTag();
    
    boost::clear_vertex( leafs_[ ch ], forest_ );
    boost::remove_vertex( leafs_[ ch ], forest_ );
    leafs_.erase( ch );
    // mother is special case ( if not isNull )
    anc_vertex j, i;
    std::map< AgentTag, anc_vertex >::iterator mother( mothers_.find( mo ) );
    if( mother != mothers_.end() ) {
        i = mother->second;
        if( i != root_ ) {
            if( boost::in_degree( i, forest_ ) == 0 ) {
                // grandma...
                j = *( boost::adjacent_vertices( i, forest_ ).first );
                
                boost::clear_vertex( i, forest_ );
                boost::remove_vertex( i, forest_ );
                mothers_.erase( mother );
                i = j;
            } else {
                // other possibility is in_degree == 1
                // if other branch is not a leaf, mother is no more a mother 
                j = *( boost::inv_adjacent_vertices( i, forest_ ).first );
                if( boost::in_degree( j, forest_ ) != 0 ) {
                    mothers_.erase( mother );
                }
            }
        }  
    } else {
        i = root_;
    }

    // further down the lineage
    while( i != root_ and boost::in_degree( i, forest_ ) == 0 ) {
        j = *( boost::adjacent_vertices( i, forest_ ).first );
        boost::clear_vertex( i, forest_ );
        boost::remove_vertex( i, forest_ );
        i = j;
    }
    // did we hit a root?
    if( i == root_ ) {
        --nr_roots_;
    }
}

void
dorsalfin::LogCsvAncestors::finalize() {
    writeForest();
}

void 
dorsalfin::LogCsvAncestors::writeHeader() {
    *log_ << "# child ( birth, x, y ) -> parent ( birth, x, y )\n";
}

void 
dorsalfin::LogCsvAncestors::writeForest() {
#ifdef DEBUG
    cout << "### Flushing!" << endl;
#endif
    // log part of the total ancestor tree, by time;
    std::map< AgentTag, anc_vertex > time_map;
    mapToTime( time_map );
    // according to map, write to log
    for( std::map< AgentTag, anc_vertex >::iterator i = time_map.begin();
        i != time_map.end(); ++i ) {
        // not printing ultimate root
        if( not i->first.isNull() ) {
            AgentTag mother =  boost::get( boost::vertex_bundle, forest_, 
                *( boost::adjacent_vertices( i->second, forest_ ).first ) );
            if( not mother.isNull() ) 
                *log_ << i->first << "\t-> " << mother << "\n";
        }
    }
}

void
dorsalfin::LogCsvAncestors::renewForest() {
    // pre: one root?
    // empty graph
    forest_.clear();
    mothers_.clear();
    vertex_nr_ = 0;
    edge_nr_ = 0;
    // add null object
    initForest();
    // and leafs at the start..
    for( std::map< AgentTag, anc_vertex >::iterator i = leafs_.begin();
        i != leafs_.end(); ++i ) {
        i->second = boost::add_vertex( forest_ );
        boost::put( boost::vertex_index, forest_, i->second, vertex_nr_++ );
        boost::put( boost::vertex_bundle, forest_, i->second, i->first );
        anc_edge aux1 = boost::add_edge( i->second, root_, forest_ ).first;
        boost::put( boost::edge_index, forest_, aux1, edge_nr_++ );
    }
    nr_roots_ = leafs_.size();
}

void
dorsalfin::LogCsvAncestors::initForest() {
    // init forest by adding the ultimate root to the mothers section
    AgentTag aux;
    root_ = boost::add_vertex( forest_ );
    boost::put( boost::vertex_index, forest_, root_, vertex_nr_++ );
    boost::put( boost::vertex_bundle, forest_, root_, aux );
    mothers_.insert( std::make_pair( aux, root_ ) );    
}

void
dorsalfin::LogCsvAncestors::mapToTime( std::map< AgentTag, anc_vertex > &tm ) {
    anc_birth_map tag_map = boost::get( boost::vertex_bundle, forest_ );
    ver_iter n, m;
    for( boost::tie( n, m ) = boost::vertices( forest_ ); n != m; ++n ) {
        tm.insert( std::make_pair( tag_map[ *n ], *n ) );
    }
}

void
dorsalfin::LogCsvAncestors::maxForestSize( uint mfs ) 
{ max_forest_size_ = mfs; }

//
// After one run, we can trace back agents to the beginning and then follow
// their development
//
dorsalfin::LogXmlAgentTrace::LogXmlAgentTrace( 
    std::string dname, std::string tracefile ) : AsyncLogObserver() {
    first_ = true;
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
    trace_ = StreamManager::instance()->openInFileStream( 
        tracefile, std::fstream::in );
    initNextTarget();
    nextTargets();
}

void 
dorsalfin::LogXmlAgentTrace::initialize( Subject *s ) { 
    update( s );
}

void
dorsalfin::LogXmlAgentTrace::update( Subject *s ) {
    GenRegAgent *ag = dynamic_cast< GenRegAgent * >( s );
    AgentTag tag = ag->myTag();
    std::set< AgentTag >::iterator hastag( targets_.find( tag ) ); 
    if( hastag != targets_.end() ) {
        // devise name
        std::string fname = tag.str();
        // write genome
        openLog( dname_ + "/" + fname + ".xml" );
        writeHeader();
#ifdef DOUBLESTRANDBREAKS
        std::vector< uint > aux( ag->nrMutations() );
        *log_ << "<mutations dsb=\"" << aux[ Chromosome::DSB ]
            << "\" cp_g=\"" << aux[ Chromosome::CP_G ] 
            << "\" rm_g=\"" << aux[ Chromosome::RM_G ] << "\"/>\n";
#endif
        *log_ << *ag;
        writeFooter();
        closeLog();
        // remove from targets_
        targets_.erase( hastag );
        // and move on
        if( targets_.empty() ) nextTargets();
    }
}

void
dorsalfin::LogXmlAgentTrace::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation dorsal_version=\"" << VERSION << "\">\n";
}

void
dorsalfin::LogXmlAgentTrace::writeFooter() {
    *log_ << "</simulation>\n";
}

void
dorsalfin::LogXmlAgentTrace::nextTargets() {
    // read all the agents with the same time of birth
    /*targets_.clear();*/
    long aux = next_target_.time;
    while( aux == next_target_.time and trace_->good() ) {
        targets_.insert( next_target_ );
        *trace_ >> next_target_.time 
                >> next_target_.loc.first 
                >> next_target_.loc.second 
                >> next_target_.i;
#ifdef DEBUG
        std::cout << next_target_ << std::endl;
#endif
    }
}

void
dorsalfin::LogXmlAgentTrace::initNextTarget() {
    // read one agent
    *trace_ >> next_target_.time >> next_target_.loc.first
            >> next_target_.loc.second >> next_target_.i;
    //std::cout << next_target_ << std::endl;
}

