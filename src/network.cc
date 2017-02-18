//
// Implementation of hopfield graph
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "network.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "gene.hh"

int dorsalfin::Network::init_propagate_ = 0;
int dorsalfin::Network::max_propagate_ = 0;
int dorsalfin::Network::nr_nodes_ = 0;

std::set< int > dorsalfin::Network::ref_ids_;
std::set< int > dorsalfin::Network::ref_input_ids_;
std::set< int > dorsalfin::Network::ref_output_ids_;
std::set< int > dorsalfin::Network::ref_transfac_ids_;

std::map< int, dorsalfin::RefNode > dorsalfin::Network::ref_nodes_;
std::map< dorsalfin::Network::dorsalfin_edge, int >
    dorsalfin::Network::ref_edges_;


dorsalfin::Network::Network() 
    : graph_(), state_( nr_nodes_, 0ul ), graph_gmap_(), graph_imap_(), 
      has_output_( false ), in_attractor_( false ) {}

dorsalfin::Network::Network( const Network &hg ) 
    : graph_(), state_( nr_nodes_, 0ul ), graph_gmap_(), graph_imap_() {
    copy( hg );
}

dorsalfin::Network::~Network() {}

dorsalfin::Network* 
dorsalfin::Network::clone() const {
    return new Network( *this );
}

void 
dorsalfin::Network::copy( const Network &hg ) {
    // not copying the graph or the two mappings, as these need to be rebuilt
    state_ = hg.state_;
    output_path_ = hg.output_path_;
    has_output_ = hg.has_output_;
    in_attractor_ = hg.in_attractor_;
}

void
dorsalfin::Network::build( const Genome &g ) {
    // pre: genes and interactions are initialised
    typedef std::multimap< int, dorsalfin_graph::vertex_descriptor > vertex_map;
    typedef std::map< Gene*, dorsalfin_graph::vertex_descriptor > gene_map;
    
    vertex_map tag_node_map;
    gene_map gene_node_map;
    std::set< int > gene_set;
    // get all genes by id in multimap
    graph_.clear();
    std::list< Chromosome* > ch = g.chromosomes();
    for( Genome::chr_iter ii = ch.begin(); ii != ch.end(); ++ii ) {
        // usually only one chromosome
        std::list< ChromosomeElement* > ces = ( **ii ).elements();
        for( Chromosome::ce_iter jj = ces.begin(); jj != ces.end(); ++jj ) {
            Gene *gg = dynamic_cast< Gene* >( *jj );
            if( gg ) {
                dorsalfin_graph::vertex_descriptor aux = 
                    boost::add_vertex( gg, graph_ );
                tag_node_map.insert( std::make_pair( gg->tag(), aux ) );
                gene_node_map.insert( std::make_pair( gg, aux ) );
                gene_set.insert( gg->tag() );
            }
        }
    }

    // check if we have an output gene
    if( std::includes( gene_set.begin(), gene_set.end(),
        ref_output_ids_.begin(), ref_output_ids_.end() ) ) {
        has_output_ = true;
        // each gene is preceded by its bindingsites, however now we also have
        // to let repeat elements break up the promotor region...
        // I should not use ref_edges_ here..
        for( Genome::chr_iter ii = ch.begin(); ii != ch.end(); ++ii ) {
            // usually only one chromosome
            bool adjacent = true;
            Gene* target = 0;
            std::list< ChromosomeElement* > ces = ( **ii ).elements();
            for( Chromosome::ce_riter j = ces.rbegin(); j != ces.rend(); ++j ) {
                if( Chromosome::IsGene()( *j ) ) {
                    adjacent = true;
                    target = dynamic_cast< Gene* >( *j );
                } else if( Chromosome::IsInteraction()( *j ) and adjacent ) {
                    Interaction *ia = dynamic_cast< Interaction* >( *j );
                    vertex_map::iterator a, b, k;
                    boost::tie( a, b ) = tag_node_map.equal_range( ia->tag() );
                    for( k = a; k != b; ++k ) {
                        boost::add_edge( k->second, gene_node_map[ target ], 
                            ia->weight(), graph_ );
                    }
#ifdef DOUBLESTRANDBREAKS
                } else if( Chromosome::IsRepeat()( *j ) ) {
                    // repeat element breaks promotor region
                    adjacent = false;
#endif
                }
            }
        }
        // cache the gene and interaction maps, not gonna change until rebuilt
        graph_gmap_ = boost::get( boost::vertex_name, graph_ );
        graph_imap_ = boost::get( boost::edge_weight, graph_ );
        // state to correct size
        state_.resize( ref_ids_.size() );
    } else {
        has_output_ = false;
    }
}

void
dorsalfin::Network::propagate() {
    std::map< Gene*, int > config;
    // can even cache these, not?
    // update as long as configurations change (so not keeping track of cycles,
    // but the max_propagate_ shouldn't be a large value anyway) and the output
    // matches the input
    
    // to store the output of the gene...
    int out = 0;
    // erase previous path
    output_path_.clear();
    bool correct = true;
    // do all the propagation
    int k = 0, l = max_propagate_;
    // adjust length of propagation to length of input string
#ifdef NOROTATE
    if( input_path_.size() + init_propagate_ < static_cast< uint >( l ) ) {
        l = input_path_.size() + init_propagate_;
    }
#endif
    while( k != l ) {
        // reset accumulated output
        out = 0;
        // store old configuration in \c config
        node_iter i, end;
        for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
            Gene *aux = boost::get( graph_gmap_, *i );
            config[ aux ] = aux->state();
        }
        // update nodes in graph according to state in copy
        // and check if any of them really changes (so we may stop)
        bool any_change = false;
        for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
            double sum = 0.0;
            edge_iter j, ant;
            for( boost::tie( j, ant ) = boost::in_edges( *i, graph_ ); 
                j != ant; ++j ) {
                Gene *aux = boost::get( graph_gmap_, 
                    boost::source( *j, graph_ ) );
                sum += config[ aux ] * boost::get( graph_imap_, *j );
            }
            Gene *aux = boost::get( graph_gmap_, *i );
            // boolean as int gives: if true then 1 else 0
            if( !close_to( sum, static_cast< double >( aux->threshold() ) ) ) {
                aux->state( static_cast< int >( sum > aux->threshold() ) );
                any_change |= ( config[ aux ] != aux->state() );
            }
            // capture output (we know it is only one output gene)
            if( *( ref_output_ids_.begin() ) == aux->tag() ) { 
                out += aux->state();
            }
        }
        // add output to path
        if( k > init_propagate_ - 1 ) {
            bool aux = ( out > 0.5 );
            correct = input_path_.test( k - init_propagate_ ) == aux;
            if( correct ) output_path_.push_back( aux );
        }
        // proceed or stop...
        if( any_change and correct ) {
            in_attractor_ = false;
            ++k;
        } else {
            in_attractor_ = not any_change;
            l = k;
        }
    }
    // but if in attractor, we still need to add more output
    if( in_attractor_ and correct ) {
        ++k;
        l = max_propagate_;
#ifdef NOROTATE
        if( input_path_.size() + init_propagate_ < static_cast< uint >( l ) ) {
            l = input_path_.size() + init_propagate_;
        }
#endif
        while( k != l ) {
            if( k > init_propagate_ - 1 ) {
                bool aux = ( out > 0.5 );
                if( input_path_.test( k - init_propagate_ ) == aux ) {
                    output_path_.push_back( aux );
                    ++k;
                } else {
                    l = k;
                }
            } else {
                ++k;
            }
        }
    }
    // output_path_ is as long as it correctly reproduced the input
}

int
dorsalfin::Network::hamming( const boost::dynamic_bitset<> &b ) const {
    // binary hamming distance between two network states (not configurations!)
    boost::dynamic_bitset<> aux( state_ );    
    // note: sizes should be equal.. if not, count the ones in the extra part
    // of our state (attractor state is zero in this region)
    // BUT only if they are not transcription factors
    int extra = 0;
    if( b.size() < aux.size() ) {
        for( uint i = b.size(); i != aux.size(); ++i ) {
            if( ref_transfac_ids_.find( i ) == ref_transfac_ids_.end() and
                aux.test( i ) ) ++extra;
        }
    }
    aux.resize( b.size() );
    aux ^= b;
    // set to zero if it is a transcription factor
    for( uint i = 0; i != aux.size(); ++i ) {
        if( ref_transfac_ids_.find( i ) != ref_transfac_ids_.end() ) {
            aux.reset( i );
        }
    }
    return aux.count() + extra;
}

void
dorsalfin::Network::doState() {
    // pre: vertex ids element of [0..n)
    state_.reset();
    node_iter i, end;
    for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
        Gene *aux = boost::get( graph_gmap_, *i );
        if( aux->state() > 0.5 ) state_.set( aux->tag() );
    }
}

void
dorsalfin::Network::setInputState( const boost::dynamic_bitset<> &b ) {
    node_iter i, end;
    for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
        Gene *aux = boost::get( graph_gmap_, *i );
        if( ref_input_ids_.find( aux->tag() ) != ref_input_ids_.end() ) {
#ifdef NOROTATE
            if( b.size() > static_cast< uint >( aux->tag() ) ) {
                aux->state( b.test( aux->tag() ) );
            } else {
                aux->state( 0 );
            }
#else
            aux->state( b.test( aux->tag() ) );
#endif
        }
    }    
}

void
dorsalfin::Network::resetState() {
    node_iter i, end;
    for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
        Gene *aux = boost::get( graph_gmap_, *i );
        aux->state( 0 );
    }    
}

bool
dorsalfin::Network::allZero() const {
    return state_.none();
}    

void
dorsalfin::Network::writeXml( std::ostream &os ) const {
    // NOTE: be aware that the id's of this graph do not match the reference
    // graph, it is a permutation...
    os << "<network>\n";
    boost::dynamic_properties dp;
    dp.property( "node_id", boost::get( boost::vertex_index, graph_ ) );
    dp.property( "iaction", boost::get( boost::edge_weight, graph_ ) );
    boost::write_graphviz( os, graph_, dp );
    os << "</network>\n";
}

void
dorsalfin::Network::readReferenceGraph( std::istream &is ) {
    ref_graph ref;
    boost::dynamic_properties dp;
    doReadReferenceGraph( is, ref, dp );
    
    // now create set of reference ids and nodes
    {
        ref_node_iter ii, end;
        for( boost::tie( ii, end ) = boost::vertices( ref ); ii != end; ++ii ) {
            // get node id and store it
            int aux = boost::get( &RefNode::node_id, ref, *ii );
            // if input
            if( boost::get( &RefNode::input, ref, *ii ) ) {
                ref_input_ids_.insert( aux );
            }
            // if output
            if( not boost::get( &RefNode::transfac, ref, *ii ) ) {
                ref_output_ids_.insert( aux );
            } else {
                // is transfac at least
                ref_transfac_ids_.insert( aux );
            }
            // and add to aggregate anyway
            ref_ids_.insert( aux );
            ref_nodes_.insert( std::make_pair( 
                aux, boost::get( boost::vertex_bundle, ref, *ii ) ) );
        }
    }
    
    // how many nodes?
    nr_nodes_ = ref_ids_.size();
    // and set of reference edges
    {
        ref_edge_iter ii, end;
        for( boost::tie( ii, end ) = boost::edges( ref ); ii != end; ++ii ) {
            // ref_edge is mapping from edge to interaction
            ref_edges_[ std::make_pair(
                boost::source( *ii, ref ), boost::target( *ii, ref ) ) ] = 
                boost::get( boost::edge_weight, ref, *ii );
        }
    }
}

void
dorsalfin::Network::initPropagate( int t ) {
    init_propagate_ = t;
}

void
dorsalfin::Network::maxPropagate( int t ) {
    max_propagate_ = t;
}

int
dorsalfin::Network::maxPropagate() {
    return max_propagate_;
}

int
dorsalfin::Network::referenceNrInputIds() {
    return ref_input_ids_.size();
}

const std::set< int > &
dorsalfin::Network::referenceIds() {
    return ref_ids_;
}

const std::set< int > &
dorsalfin::Network::referenceInputIds() {
    return ref_input_ids_;
}

const std::set< int > &
dorsalfin::Network::referenceTransFacIds() {
    return ref_transfac_ids_;
}

const std::map< int, dorsalfin::RefNode > &
dorsalfin::Network::referenceNodes() {
    return ref_nodes_;
}

const std::map< dorsalfin::Network::dorsalfin_edge, int > &
dorsalfin::Network::referenceEdges() {
    return ref_edges_;
}

void
dorsalfin::Network::writeReferenceGraph( std::ostream &os ) {
    typedef std::map< int, ref_graph::vertex_descriptor > vertex_map;
    
    ref_graph ref;
    vertex_map gene_node_map;
    // build graph from ref_nodes_ and ref_edges_
    for( std::map< int, RefNode >::const_iterator ii = ref_nodes_.begin();
        ii != ref_nodes_.end(); ++ii ) {
        gene_node_map.insert( 
            std::make_pair( ii->first, boost::add_vertex( ii->second, ref ) ) );
    }
    // and the edges
    for( edge_iaction_citer i = ref_edges_.begin(); i != ref_edges_.end(); ++i){
        boost::add_edge( gene_node_map[ i->first.first ],
            gene_node_map[ i->first.second ], i->second, ref );
    }

    boost::dynamic_properties dp;
    dp.property( "node_id", boost::get( &RefNode::node_id, ref ) );
    dp.property( "thr", boost::get( &RefNode::threshold, ref ) );
    dp.property( "transfac", boost::get( &RefNode::transfac, ref ) );
    dp.property( "input", boost::get( &RefNode::input, ref ) );
    dp.property( "iaction", boost::get( boost::edge_weight, ref ) );
    boost::write_graphviz( os, ref, dp );
}    

void
dorsalfin::Network::doReadReferenceGraph( std::istream &is, ref_graph &ref,
    boost::dynamic_properties &dp ) {
    // auxilary function
    dp.property( "node_id", boost::get( &RefNode::node_id, ref ) );
    dp.property( "thr", boost::get( &RefNode::threshold, ref ) );
    dp.property( "transfac", boost::get( &RefNode::transfac, ref ) );
    dp.property( "input", boost::get( &RefNode::input, ref ) );
    dp.property( "iaction", boost::get( boost::edge_weight, ref ) );

    try {
        boost::read_graphviz( is, ref, dp );
    } catch( const boost::graph_exception &e ) {
        std::cout << "Exception: " << e.what() << "\n";
    } catch( const boost::dynamic_property_exception &e ) {
        std::cout << "Exception: " << e.what() << "\n";
    }
}


