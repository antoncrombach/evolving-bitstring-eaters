//
// Chromosome implementation...
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "pool.hh"
#include "ref_node.hh"
#include "network.hh"
#include "chromosome.hh"


template<> dorsalfin::ObjectCache< dorsalfin::Chromosome >* 
dorsalfin::ObjectCache< dorsalfin::Chromosome >::instance_ = 0;

double dorsalfin::Chromosome::cp_tp_rate_ = 0.0;
double dorsalfin::Chromosome::rm_tp_rate_ = 0.0;
double dorsalfin::Chromosome::rm_ltr_rate_ = 0.0;
double dorsalfin::Chromosome::new_tp_rate_ = 0.0;
double dorsalfin::Chromosome::dsb_recombination_ = 0.0;

double dorsalfin::Chromosome::cp_gene_rate_ = 0.0;
double dorsalfin::Chromosome::rm_gene_rate_ = 0.0;
double dorsalfin::Chromosome::thr_gene_rate_ = 0.0;
double dorsalfin::Chromosome::tag_gene_rate_ = 0.0;
double dorsalfin::Chromosome::state_rate_ = 0.0;

double dorsalfin::Chromosome::cp_ia_rate_ = 0.0;
double dorsalfin::Chromosome::rm_ia_rate_ = 0.0;
double dorsalfin::Chromosome::new_ia_rate_ = 0.0;
double dorsalfin::Chromosome::tag_ia_rate_ = 0.0;
double dorsalfin::Chromosome::weight_ia_rate_ = 0.0;

// constant
int dorsalfin::Chromosome::NR_MUT = 15;

// note: using magic number
dorsalfin::Chromosome::Chromosome() 
    : parent_( 0 ), chro_( new std::list< ChromosomeElement* >() ),
      mut_events_( NR_MUT, 0 ), retros_( 0 ), ltrs_( 0 ), len_( 0 ),
      update_retro_( true ), update_ltr_( true ), update_len_( true ) {}

// note: using magic number
dorsalfin::Chromosome::Chromosome( Genome *g, 
    std::list< ChromosomeElement* > *ll ) 
    : parent_( g ), chro_( ll ), mut_events_( NR_MUT, 0 ), retros_( 0 ),
      ltrs_( 0 ), len_( ll->size() ), update_retro_( true ),
      update_ltr_( true ), update_len_( true ) {}

dorsalfin::Chromosome::Chromosome( const Chromosome &c ) 
    : chro_( new std::list< ChromosomeElement* >() ), mut_events_( NR_MUT, 0 ) {
    copy( c );
}

dorsalfin::Chromosome::~Chromosome() {
    smart_erase( *chro_, chro_->begin(), chro_->end() );
    delete chro_;
}

void
dorsalfin::Chromosome::copy( const Chromosome &c ) {
    parent_ = c.parent_;
    for( ce_iter i = c.chro_->begin(); i != c.chro_->end(); ++i ) {
        ChromosomeElement *aux = ( **i ).clone();
        chro_->push_back( aux );
    }
    std::copy( c.mut_events_.begin(), c.mut_events_.end(), mut_events_.begin());
    retros_ = c.retros_;
    ltrs_ = c.ltrs_;
    len_ = c.len_;
    update_retro_ = c.update_retro_;
    update_ltr_ = c.update_ltr_;
    update_len_ = c.update_len_;
}

dorsalfin::Chromosome*
dorsalfin::Chromosome::clone() const {
    Chromosome *cc = ObjectCache< Chromosome >::instance()->borrowObject();
    cc->copy( *this );
    return cc;
}

bool
dorsalfin::Chromosome::toPool() {
    // empty everything
    smart_return( *chro_, chro_->begin(), chro_->end() );
    std::fill_n( mut_events_.begin(), NR_MUT, 0 );
    retros_ = 0;
    ltrs_ = 0;
    len_ = 0;    
    update_retro_ = true;
    update_ltr_ = true;
    update_len_ = true;
    return ObjectCache< Chromosome >::instance()->returnObject( this );
}

//
//
// Adding some interactions, building up the genome/chromosome from a network
//
//
void
dorsalfin::Chromosome::initialise() {
    // pre: assuming genes have ids from 0 to n in the dot-file
    std::map< int, RefNode > ref_nodes = Network::referenceNodes();
    std::map< Network::dorsalfin_edge, int > ref_edges =
        Network::referenceEdges();
    std::multimap< int, int > reverse_edges;
    for( std::map< Network::dorsalfin_edge, int >::iterator i =
        ref_edges.begin(); i != ref_edges.end(); ++i ) {
        reverse_edges.insert( reverse( i->first ) );
    }
    
    for( ce_iter jj = chro_->begin(); jj != chro_->end(); ++jj ) {
        Gene *gg = dynamic_cast< Gene* >( *jj );
        if( gg ) {
            gg->import( ref_nodes[ gg->tag() ] );
            std::multimap< int, int >::const_iterator i, j, k;
            boost::tie( i, j ) = reverse_edges.equal_range( gg->tag() );
            for( k = i; k != j; ++k ) {
                // what about parallel edges? Not allowed here!
                chro_->insert( jj, new Interaction( k->second, 
                    ref_edges.find( reverse( *k ) )->second ) );
            }
        }
    }
    update_len_ = true;
}

//
//
// Mutation related methods
//
//
int
dorsalfin::Chromosome::mutate() {
    // pre: all chromosomes have been reset before mutating them all
    // iterate along chromosomes
    ce_iter i = chro_->begin();
    while( i != chro_->end() ) {
        // note: it is a bit tricky, but incrementing i is done in the methods
        if( ( **i ).isActive() ) {
            if( IsInteraction()( *i ) ) {
                i = interactionMutate( i );
            } else if( IsGene()( *i ) ) {
                i = geneMutate( i );
            }
#ifdef DOUBLESTRANDBREAKS
              else if( IsRepeat()( *i ) ) {
                i = repeatMutate( i );
            } else if( IsRetroposon()( *i ) ) {
                i = retroposonMutate( i );
            }
#endif
        } else {
            ++i;
        }
    }
    // and maybe something new arises
#ifdef DOUBLESTRANDBREAKS
    retroposonInvade();
#endif
    interactionGenesis();
    return 0;
}

void
dorsalfin::Chromosome::retroposonInvade() {
    if( close_to( new_tp_rate_, 0.0 ) ) {
        return;
    }
    // assuming small rate
    int hux = size();
    if( uniform() < 1.0 - pow( 1.0 - new_tp_rate_, hux ) ) {
        // find a spot
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeLocation();
        // insert
        std::list< ChromosomeElement* > eux;
        ChromosomeElement *dux = new Repeat();
        dux->inactivate();
        eux.push_back( dux );
        // FIX using magic number!!!   
        dux = new Retroposon( 1000 );
        dux->inactivate();
        eux.push_back( dux );
        dux = new Repeat();
        dux->inactivate();
        eux.push_back( dux );
        aux->splice( bux, eux, 3, 1 );
    }
}

void
dorsalfin::Chromosome::interactionGenesis() {
    if( close_to( new_ia_rate_, 0.0 ) ) {
        return;
    }
    // assuming small rate
    int hux = size();
    if( uniform() < 1.0 - pow( 1.0 - new_ia_rate_, hux ) ) {
        // find a spot (always inserts in front of a random gene)
        Chromosome *aux;
        ce_iter bux, cux;
        boost::tie( aux, bux ) = randGenomeGene();
        // get tag and weight
        std::set< int > ref_ids = Network::referenceTransFacIds();
        int t = *( random_element( ref_ids.begin(), ref_ids.end(), 
            rand_range< int > ) );
        int w = 2 * static_cast< int >( uniform() + 0.5 ) - 1;
        // insert
        ChromosomeElement *dux = new Interaction( t, w );
        dux->inactivate();
        aux->insert( bux, dux );
        ++aux->mut_events_[ NW_IA ];
    }
}

dorsalfin::Chromosome::ce_iter
dorsalfin::Chromosome::retroposonMutate( ce_iter i ) {
    // pre: i is active and retroposon
    if( close_to( cp_tp_rate_ + rm_tp_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    double uu = uniform();
    if( uu < cp_tp_rate_ ) {
        // insert a copy of the retroposon and LTRs somewhere
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeLocation();
        // get the repeats
        ce_iter ll( boost::prior( i ) );
        ce_iter rr( boost::next( i, 2 ) );
        std::list< ChromosomeElement* > cux;
        for( ce_iter j = ll; j != rr; ++j ) {
            ChromosomeElement *dux = ( **j ).clone();
            dux->inactivate();
            cux.push_back( dux );
        }
        aux->splice( bux, cux, 3, 1 );
        ++aux->mut_events_[ CP_RP ];
        ++i;
    } else if( uu < cp_tp_rate_ + rm_tp_rate_ ) {
        // reciprocal recombination, one LTR stays
        ce_iter ll( boost::prior( i ) );
        ce_iter rr( boost::next( i ) );
        i = smart_return( *chro_, ll, rr );
        // inactivate leftover repeat
        ( **i ).inactivate();
        ++mut_events_[ RM_RP ];
        --retros_;
        --ltrs_;
        len_ -= 2;
    } else {
        ++i;
    }
    return i;
}

dorsalfin::Chromosome::ce_iter
dorsalfin::Chromosome::repeatMutate( ce_iter i ) {
    if( close_to( dsb_recombination_ + rm_ltr_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    Repeat *aux = dynamic_cast< Repeat* >( *i );
    double uu = uniform();
    if( uu < dsb_recombination_ ) {
        // needs to be repaired @ genome level
        aux->induceDSB();
        ( **i ).inactivate();
        ++i;
        ++mut_events_[ DSB ];
    } else if( uu < dsb_recombination_ + rm_ltr_rate_ ) {
        // first check if flanking a retroposon
        if( len_ > 1 ) {
            if( i == chro_->begin() ) {
                if( IsRetroposon()( *( boost::next( i ) ) ) ) {
                    ++i;
                } else {
                    i = smart_return( *chro_, i );
                    ++mut_events_[ RM_LTR ];
                    --ltrs_;
                    --len_;
                }
            } else if( boost::next( i ) == chro_->end() ) {
                if( IsRetroposon()( *( boost::prior( i ) ) ) ) {
                    ++i;
                } else {
                    i = smart_return( *chro_, i );
                    ++mut_events_[ RM_LTR ];
                    --ltrs_;
                    --len_;
                }
            } else {
                if( IsRetroposon()( *( boost::prior( i ) ) ) or
                    IsRetroposon()( *( boost::next( i ) ) ) ) {
                    ++i;
                } else {
                    i = smart_return( *chro_, i );
                    ++mut_events_[ RM_LTR ];
                    --ltrs_;
                    --len_;
                }
            }
        } else {
            i = smart_return( *chro_, i );
            ++mut_events_[ RM_LTR ];
            --ltrs_;
            --len_;
        }            
    } else {
        // just in case (patch)
        aux->repairDSB();
        ++i;
    }
    return i;
}


dorsalfin::Chromosome::ce_iter
dorsalfin::Chromosome::geneMutate( ce_iter i ) {
    // pre: i is active
    if( close_to( cp_gene_rate_ + rm_gene_rate_ + thr_gene_rate_ +
        tag_gene_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    double uu = uniform();
    if( uu < cp_gene_rate_ ) {
        // insert a copy of the current gene somewhere in the genome
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeGene();
        // find the gene
        ce_iter ll( interactionSelect( i ) );
        ce_iter rr( boost::next( i ) );
        std::list< ChromosomeElement* > cux;
        uint tt = 0;
        for( ce_iter j = ll; j != rr; ++j ) {
            ChromosomeElement *dux = ( **j ).clone();
            dux->inactivate();
            cux.push_back( dux );
            ++tt;
        }
        aux->splice( boost::next( bux ), cux, tt, 0 );
        ++aux->mut_events_[ CP_G ];
        ++i;
    } else if( uu < cp_gene_rate_ + rm_gene_rate_ ) {
        // delete the current gene and its interactions
        ce_iter ll( interactionSelect( i ) );
        ce_iter rr( boost::next( i ) );
        len_ -= std::distance( ll, rr );
        i = smart_return( *chro_, ll, rr );
        ++mut_events_[ RM_G ];
    } else if( uu < cp_gene_rate_ + rm_gene_rate_ + thr_gene_rate_ ) {
        dynamic_cast< Gene* >( *i )->mutateThreshold();
        ++mut_events_[ THR ];
        ++i;
    } else if( uu < cp_gene_rate_ + rm_gene_rate_ + thr_gene_rate_ +
        tag_gene_rate_ ) {
        // change target to some value.. pick a random gene :)
        std::set< int > ref = Network::referenceIds();
        int t = *( random_element( ref.begin(), ref.end(), 
            rand_range< int > ) );
        dynamic_cast< Gene* >( *i )->tag( t );
        ++mut_events_[ TAG ];
        ++i;
    } else {
        ++i;
    }
    return i;
}

dorsalfin::Chromosome::ce_iter
dorsalfin::Chromosome::interactionMutate( ce_iter i ) {
    if( close_to( cp_ia_rate_ + rm_ia_rate_ + 
        tag_ia_rate_ + weight_ia_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    double uu = uniform();
    if( uu < weight_ia_rate_ ) {
        dynamic_cast< Interaction* >( *i )->toggle();
        ++mut_events_[ W_IA ];
        ++i;
    } else if( uu < weight_ia_rate_ + cp_ia_rate_ ) {
        // copy this interaction
        Chromosome *bux;
        ce_iter cux;
        boost::tie( bux, cux ) = randGenomeGene();
        ce_iter j( bux->insert( cux, ( **i ).clone() ) );
        ( **j ).inactivate();
        ++bux->mut_events_[ CP_IA ];
        ++i;
    } else if( uu < weight_ia_rate_ + cp_ia_rate_ + rm_ia_rate_ ) {
        // remove this interaction
        i = smart_return( *chro_, i );
        ++mut_events_[ RM_IA ];
        --len_;
    } else if( uu < weight_ia_rate_ + cp_ia_rate_ + 
        rm_ia_rate_ + tag_ia_rate_ ) {
        // change target to some value.. pick a random gene :)
        std::set< int > ref = Network::referenceTransFacIds();
        int t = *( random_element( ref.begin(), ref.end(), 
            rand_range< int > ) );
        dynamic_cast< Interaction* >( *i )->tag( t );
        ++i;
        ++mut_events_[ T_IA ];
    } else {
        ++i;
    }
    return i;
}

int
dorsalfin::Chromosome::genePerturb() {
    if( close_to( state_rate_, 0.0 ) ) {
        return 0;
    }
    // assuming small rate
    int aux = size(), fux = 0;
    if( uniform() < 1.0 - pow( 1.0 - state_rate_, aux ) ) {
        // find a spot
        Chromosome *bux;
        ce_iter cux;
        boost::tie( bux, cux ) = randGenomeGene();
        dynamic_cast< Gene* >( *cux )->mutateState();
        ++mut_events_[ ST ];
        ++fux;
    }
    return fux;
}

dorsalfin::Chromosome::ce_iter
dorsalfin::Chromosome::interactionSelect( ce_iter i ) {
    // pre: IsGene( i )
    ce_riter ii( i );
    ce_riter jj = chro_->rend();
    // BLS for first non-binding site
    while( ii != jj ) {
        if( IsInteraction()( *ii ) ) {
            ++ii;
        } else {
            jj = ii;
        }
    }
    return jj.base();
}

//
//
// Genome related methods
//
//
std::list< dorsalfin::Chromosome* >
dorsalfin::Chromosome::segments() {
    std::list< Chromosome* > result;
    if( mut_events_[ DSB ] == 0 ) {
        // creating an alias...if not empty list
        if( !chro_->empty() ) {
            result.push_back( this );
        } else {
            Chromosome *cux = 
                ObjectCache< Chromosome >::instance()->borrowObject();
            cux->parent_ = parent_;
            result.push_back( cux );
        }
    } else {
        // creation of new chromosomes automatically sets their update flag
        // get the bits and pieces
        bool found = false;
        ce_iter jj( chro_->begin() );
        ce_iter ii( jj );
        // and loop...
        while( jj != chro_->end() ) {
            Repeat *aux = dynamic_cast< Repeat* >( *jj );
            if( aux ) {
                if( aux->hasDSB() ) {
                    aux->repairDSB();
                    found = true;
                }
            }
            ++jj;
            // if found, put the segment in the vector
            if( found ) {
                Chromosome *cux = 
                    ObjectCache< Chromosome >::instance()->borrowObject();
                cux->parent_ = parent_;
                cux->chro_->splice( cux->chro_->end(), *chro_, ii, jj );
                result.push_back( cux );
                ii = jj;
                found = false;
            }
        }
        // and the last one
        Chromosome *cux = 
            ObjectCache< Chromosome >::instance()->borrowObject();
        cux->parent_ = parent_;
        cux->chro_->splice( cux->chro_->end(), *chro_, ii, jj );
        result.push_back( cux );
    }
    // an empty chromosome is the leftover, if there were segments
    return result;
}

void
dorsalfin::Chromosome::append( Chromosome *chr ) {
    // not copying parent_
    chro_->splice( chro_->end(), *( chr->chro_ ) );
}

void
dorsalfin::Chromosome::reset() {
    // pre: chro_ is initialised
    for( ce_iter i = chro_->begin(); i != chro_->end(); ++i ) {
        ( **i ).activate();
    }
    // making sure all entries exist and are zero
    std::fill_n( mut_events_.begin(), NR_MUT, 0 );
}

void
dorsalfin::Chromosome::recache() {
    update_retro_ = true;
    update_ltr_ = true;
    update_len_ = true;
}

//
//
// Private auxilary functions
//
//
dorsalfin::Chromosome::ce_iter
dorsalfin::Chromosome::insert( ce_iter i, ChromosomeElement *ce ) {
    ++len_;
    return chro_->insert( i, ce );
}

void
dorsalfin::Chromosome::splice( ce_iter i, 
    std::list< ChromosomeElement* > ces, uint ll, uint rr ) {
    // overloading list method to update cached length and retroposons
    len_ += ll;
    retros_ += rr;
    ltrs_ += 2 * rr;
    chro_->splice( i, ces );
}

//
//
// Random location/gene
//
//
boost::tuple< dorsalfin::Chromosome*, dorsalfin::Chromosome::ce_iter >
dorsalfin::Chromosome::randChromosomeLocation() {
    // pre: chro_ is initialised
    std::list< ce_iter > aux;
    ce_iter i = chro_->begin();
    if( IsGene()( *i ) or IsRepeat()( *i ) ) {
        aux.push_back( i );
    }
    ++i;
    while( i != chro_->end() ) {
        if( IsGene()( *i ) ) {
            aux.push_back( i );
        } else if( IsRepeat()( *i ) and 
                   !IsRetroposon()( *( boost::prior( i ) ) ) ) {
            aux.push_back( i );
        }
        ++i;
    }
    // shouldn't I be adding chro_->end() Yes :)
    aux.push_back( chro_->end() );
    return boost::make_tuple( this, 
        *( random_element( aux.begin(), aux.end(), rand_range< int > ) ) );
}

boost::tuple< dorsalfin::Chromosome*, dorsalfin::Chromosome::ce_iter >
dorsalfin::Chromosome::randChromosomeGene() {
    // pre: chro_ is initialised
    std::list< ce_iter > aux;
    ce_iter i = chro_->begin();
    while( i != chro_->end() ) {
        if( IsGene()( *i ) ) {
            aux.push_back( i );
        }
        ++i;
    }
    // if there is no gene, return chro_->end()
    if( not aux.empty() ) {
        return boost::make_tuple( this, 
            *( random_element( aux.begin(), aux.end(), rand_range< int > ) ) );
    } else {
        return boost::make_tuple( this, chro_->end() );
    }
}

//
//
// Set the state of the network
//
//
void
dorsalfin::Chromosome::setState( const boost::dynamic_bitset<> &state ) {
    for( ce_iter i = chro_->begin(); i != chro_->end(); ++i ) {
        Gene *aux = dynamic_cast< Gene* >( *i );
        if( aux ) {
            if( aux->tag() < static_cast< int >( state.size() ) ) {
                aux->state( static_cast< int >( state[ aux->tag() ] ) );
            } else {
                aux->state( 0 );
            }
        }
    }
}

void
dorsalfin::Chromosome::setState( int tag, int state ) {
    for( ce_iter i = chro_->begin(); i != chro_->end(); ++i ) {
        Gene *aux = dynamic_cast< Gene* >( *i );
        if( aux ) {
            if( aux->tag() == tag ) {
                aux->state( state );
            }
        }
    }
}

//
//
// Get/Set some properties of chromosomes
//
//
uint
dorsalfin::Chromosome::nrRepeats() const {
    if( update_ltr_ ) {
        update_ltr_ = false;
        ltrs_ = std::count_if( chro_->begin(), chro_->end(), IsRepeat() ); 
    }
    return ltrs_;
}

uint
dorsalfin::Chromosome::nrRetroposons() const {
    if( update_retro_ ) {
        update_retro_ = false;
        retros_ = std::count_if( chro_->begin(), chro_->end(),
            IsRetroposon() );
    }
    return retros_;
}

uint
dorsalfin::Chromosome::size() const { 
    if( update_len_ ) {
        update_len_ = false;
        len_ = chro_->size();
    }
    return len_;
}

//
//
// Writing the chromosome to a stream
//
//
void
dorsalfin::Chromosome::writeXml( std::ostream &os ) const {
    // pre: chro_ is initialised
    os << "<chromosome len=\"" << len_ << "\">\n";
    ce_iter i = chro_->begin(); 
    while( i != chro_->end() ) {
        os << **i;
        ++i;
    }
    os << "</chromosome>\n";
}

//
//
// Zero the mutation rates (static method)
//
//
void
dorsalfin::Chromosome::zeroRates() {
    cp_tp_rate_ = 0.0;
    rm_tp_rate_ = 0.0;
    rm_ltr_rate_ = 0.0;
    new_tp_rate_ = 0.0;
    dsb_recombination_ = 0.0;
    
    cp_gene_rate_ = 0.0;
    rm_gene_rate_ = 0.0;
    thr_gene_rate_ = 0.0;
    tag_gene_rate_ = 0.0;
    state_rate_ = 0.0;
    
    cp_ia_rate_ = 0.0;
    rm_ia_rate_ = 0.0;
    new_ia_rate_ = 0.0;
    tag_ia_rate_ = 0.0;
    weight_ia_rate_ = 0.0;
}

//
//
// Static methods for all the parameter rates
//
//
void
dorsalfin::Chromosome::copyRetroposonRate( double f )
{ cp_tp_rate_ = f; }

double
dorsalfin::Chromosome::copyRetroposonRate() 
{ return cp_tp_rate_; }

void
dorsalfin::Chromosome::removeRetroposonRate( double f )
{ rm_tp_rate_ = f; }

double
dorsalfin::Chromosome::removeRetroposonRate()
{ return rm_tp_rate_; }

void
dorsalfin::Chromosome::removeRepeatRate( double f )
{ rm_ltr_rate_ = f; }

double
dorsalfin::Chromosome::removeRepeatRate()
{ return rm_ltr_rate_; }

void
dorsalfin::Chromosome::newRetroposonRate( double f )
{ new_tp_rate_ = f; }

double
dorsalfin::Chromosome::newRetroposonRate()
{ return new_tp_rate_; }

void 
dorsalfin::Chromosome::recombinationRate( double f )
{ dsb_recombination_ = f; }

double
dorsalfin::Chromosome::recombinationRate()
{ return dsb_recombination_; }

void
dorsalfin::Chromosome::copyGeneRate( double f ) 
{ cp_gene_rate_ = f; }

double 
dorsalfin::Chromosome::copyGeneRate()
{ return cp_gene_rate_; }

void 
dorsalfin::Chromosome::removeGeneRate( double f )
{ rm_gene_rate_ = f; }

double 
dorsalfin::Chromosome::removeGeneRate()
{ return rm_gene_rate_; }

void
dorsalfin::Chromosome::tagGeneRate( double d ) 
{ tag_gene_rate_ = d; }

void
dorsalfin::Chromosome::thresholdGeneRate( double d ) 
{ thr_gene_rate_ = d; }

void 
dorsalfin::Chromosome::weightInteractionRate( double d ) 
{ weight_ia_rate_ = d; }

void 
dorsalfin::Chromosome::copyInteractionRate( double d ) 
{ cp_ia_rate_ = d; }

void 
dorsalfin::Chromosome::removeInteractionRate( double d ) 
{ rm_ia_rate_ = d; }

void 
dorsalfin::Chromosome::newInteractionRate( double d ) 
{ new_ia_rate_ = d; }

void 
dorsalfin::Chromosome::tagInteractionRate( double d ) 
{ tag_ia_rate_ = d; }

void
dorsalfin::Chromosome::stateRate( double d )
{ state_rate_ = d; }

