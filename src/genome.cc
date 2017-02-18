//
// Genome class implementation.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "network.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "chromelement.hh"

dorsalfin::Genome::Genome() : chromos_( new std::list< Chromosome* >() ),
    nr_mutations_( Chromosome::NR_MUT, 0 ), nr_perturbed_( 0 ) {}

// note: using magic number
dorsalfin::Genome::Genome( std::list< Chromosome* > *c ) 
    : chromos_( c ), nr_mutations_( Chromosome::NR_MUT, 0 ),
      nr_perturbed_( 0 ) {
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).parent( this );
    }
}

dorsalfin::Genome::Genome( const Genome &g ) 
    : nr_mutations_( Chromosome::NR_MUT, 0 ), nr_perturbed_( 0 ) {
    Genome::copy( g );
}

dorsalfin::Genome::~Genome() {
    clear();
    delete chromos_;
}

dorsalfin::Genome*
dorsalfin::Genome::clone() const {
    return new Genome( *this );
}

void
dorsalfin::Genome::copy( const Genome &g ) {
    chromos_ = new std::list< Chromosome* >();
    for( chr_iter i = g.chromos_->begin(); i != g.chromos_->end(); ++i ) {
        chromos_->push_back( ( **i ).clone() );
        chromos_->back()->parent( this );
    }
    std::copy( g.nr_mutations_.begin(), g.nr_mutations_.end(),
        nr_mutations_.begin() );
    nr_perturbed_ = g.nr_perturbed_;
}

void
dorsalfin::Genome::clear() {
    //smart_erase( *chromos_, chromos_->begin(), chromos_->end() );
    smart_return( *chromos_, chromos_->begin(), chromos_->end() );
}

//
//
// Building genome a bit further with help of network
//
//
void
dorsalfin::Genome::initialise() {
    // pre: reference graph has been read from file
    // for each gene add interactions
    for( chr_iter ii = chromos_->begin(); ii != chromos_->end(); ++ii ) {
        // usually only one chromosome
        ( **ii ).initialise();
    }
}

//
//
// Reproduction shizzle
//
//
void
dorsalfin::Genome::duplicate() {
    std::list< Chromosome* > aux;
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        aux.push_back( ( **i ).clone() );
    }
    chromos_->splice( chromos_->end(), aux );
}

dorsalfin::Genome*
dorsalfin::Genome::split() {
    // note: assuming only two chromosomes
    std::list< Chromosome* > *aux = new std::list< Chromosome* >();
    aux->splice( aux->begin(), *chromos_, boost::next( chromos_->begin() ) );
    Genome *result = new Genome( aux );
    return result;
}

int 
dorsalfin::Genome::mutate() {
    int result = 0;
    // first reset all counters etc
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).reset();
    }
    // and then mutate the shizzle
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).mutate();
    }

    //
    // as we do not use DSBs, we can leave the recombination code out :)
    //
    // Gather a few counts, first clear
    std::fill_n( nr_mutations_.begin(), Chromosome::NR_MUT, 0 );
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        // genes
        std::vector< uint > aux = ( **i ).nrMutations();
        // might be empty if a cached version is taken...
        if( ! aux.empty() ) {
            std::transform( nr_mutations_.begin(), nr_mutations_.end(),
                aux.begin(), nr_mutations_.begin(), std::plus< uint >() );
        }
    }

#ifdef DOUBLESTRANDBREAKS
    // get all segments
    std::list< Chromosome* > recombined;
    std::vector< Chromosome* > heads, middles, tails;
    heads.reserve( chromos_->size() );
    tails.reserve( chromos_->size() );
    chr_iter i = chromos_->begin();
    while( i != chromos_->end() ) {
        // get the segments
        std::list< Chromosome* > aux = ( **i ).segments();
        // do something depending on # segments
        if( aux.size() == 1 ) {
            // whole chromosome
            recombined.push_back( aux.front() );
        } else if( aux.size() == 2 ) {
            // head and tail
            heads.push_back( aux.front() );
            tails.push_back( aux.back() );
        } else {
            // head and tail with middle stuff too
            middles.reserve( middles.size() + aux.size() - 2 );
            heads.push_back( aux.front() );
            tails.push_back( aux.back() );
            // and all the middle parts
            std::copy( boost::next( aux.begin() ), boost::prior( aux.end() ),
                std::back_inserter( middles ) );
        }
        ++i;
    }

    // try to reuse these...
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        if( ( **i ).empty() ) {
            if( !( **i ).toPool() ) delete *i;
        }
    }
    // chromosomes are empty
    chromos_->clear();

    // randomly assign middle segments to heads
    std::random_shuffle( middles.begin(), middles.end(), rand_range< int > );
    for( std::vector< Chromosome* >::iterator i = middles.begin(); 
        i != middles.end(); ++i ) {
        Chromosome *aux = *( random_element( heads.begin(),
             heads.end(), rand_range< int > ) );
        aux->append( *i );
    }
    //smart_erase( middles, middles.begin(), middles.end() );
    smart_return( middles, middles.begin(), middles.end() );

    // randomly assign the last part of the chromosomes
    std::random_shuffle( tails.begin(), tails.end(), rand_range< int > );
    for( std::vector< Chromosome* >::iterator i = heads.begin(); 
        i != heads.end(); ++i ) {
        ( **i ).append( tails.back() );
        delete tails.back();
        tails.pop_back();
    }

    // make them recache their info
    for( std::vector< Chromosome* >::iterator i = heads.begin(); 
        i != heads.end(); ++i ) {
        ( **i ).recache();
    }

    // and insert whole chromosomes into chromos_
    chromos_->splice( chromos_->begin(), recombined );
    // add the truly recombined ones
    std::copy( heads.begin(), heads.end(), std::back_inserter( *chromos_ ) );
    heads.clear();
#endif
    return result; 
}

//
//
// During lifetime we may have some perturbations (reset counter)
//
//
void
dorsalfin::Genome::perturb() {
    nr_perturbed_ = 0;
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        nr_perturbed_ += ( **i ).genePerturb();
    }
}

void
dorsalfin::Genome::setState( const boost::dynamic_bitset<> &state ) {
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).setState( state );
    }
}

void
dorsalfin::Genome::setState( int tag, int state ) {
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).setState( tag, state );
    }
}

void
dorsalfin::Genome::setState( const std::set< int > &tags, 
    const boost::dynamic_bitset<> &state ) {
    for( std::set< int >::const_iterator t = tags.begin(); t != tags.end(); 
        ++t ) {
        setState( *t, state.test( *t ) );
    }
}

//
//
// Needed by chromosomes
//
//
boost::tuple< dorsalfin::Chromosome*, 
    std::list< dorsalfin::ChromosomeElement* >::iterator >
dorsalfin::Genome::randLocation() {
    Chromosome *bux = *( random_element( chromos_->begin(), chromos_->end(),
        rand_range< int > ) );
    return bux->randChromosomeLocation();
}

boost::tuple< dorsalfin::Chromosome*, 
    std::list< dorsalfin::ChromosomeElement* >::iterator >
dorsalfin::Genome::randGene() {
    Chromosome *bux = *( random_element( chromos_->begin(), chromos_->end(),
        rand_range< int > ) );
    return bux->randChromosomeGene();
}

//
//
// Writing xml to file
//
//
void
dorsalfin::Genome::writeXml( std::ostream &os ) const {
    os << "<genome>\n";
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        os << **i;
    }
    os << "</genome>\n";
}

//
//
// Reporting some properties of the genome
//
//
uint
dorsalfin::Genome::length() const {
    uint result = 0;
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).size();
    }
    return result;
}

uint
dorsalfin::Genome::nrRetroposons() const {
    uint result = 0;
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).nrRetroposons();
    }
    return result;
}

uint
dorsalfin::Genome::nrRepeats() const {
    uint result = 0;
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).nrRepeats();
    }
    return result;
}

uint
dorsalfin::Genome::nrDoubleStrandBreaks() const {
    uint result = 0;
    for( chr_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).nrDoubleStrandBreaks();
    }
    return result;
}

std::vector< uint >
dorsalfin::Genome::nrMutations() const {
    return nr_mutations_;
}

bool
dorsalfin::Genome::hasMutation() const {
    return std::accumulate( nr_mutations_.begin(), nr_mutations_.end(), 0 ) > 0;
}
    

