//
// Implementation of a few observers that dump their data in files.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "loggenome.hh"

//#include "environment.hh"
#include "population.hh"
#include "genreg_agent.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "retroposon.hh"
#include "gene.hh"

//#include "duo_agent.hh"

//
// Simple csv observer log a bunch of genome properties
//
dorsalfin::LogGenomeStats::LogGenomeStats( long i ) 
    : LogObserver( i ), sample_( 0 ), nr_repeats_( 0 ), nr_retroposons_( 0 ), 
        nr_genes_(), nr_bsites_(), in_degree_(), out_degree_() {}

void
dorsalfin::LogGenomeStats::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the configurations, bsites and genes
    // not completely sure if this is OOP wise, but it will parallelize better
    nr_genes_.clear();
    nr_bsites_.clear();
    in_degree_.clear();
    out_degree_.clear();
#ifdef DOUBLESTRANDBREAKS
    nr_repeats_ = 0;
    nr_retroposons_ = 0;
#endif

    sample_ = 0;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        if( na ) {
            std::list< Chromosome* > chr( na->genome()->chromosomes() );
            // should be only one chromosome here...
            /*for( Genome::chr_iter i = chr.begin(); i != chr.end(); ++i ) {
                doCounting( *i );
            }*/
            doCounting( chr.front() );
            ++sample_;
        }
    }
    // and tell the child observers...
    notify();
}

void
dorsalfin::LogGenomeStats::doCounting( Chromosome *ch ) {
    // make histograms for one chromosome, add final histo to the distribution
    histo aux_genes, aux_bsites, aux_indegree, aux_outdegree;
    std::set< int > tags;
    std::list< ChromosomeElement* > ces( ch->elements() );
    for( Chromosome::ce_iter i = ces.begin(); i != ces.end(); ++i ) {
        if( Chromosome::IsGene()( *i ) ) {
            Gene *aux = dynamic_cast< Gene* >( *i );
            // gene cp number
            aux_genes[ aux->tag() ] += 1;
            // indegree
            aux_indegree[ aux->tag() ] += tags.size();
            // outdegree
            for( std::set< int >::iterator i = tags.begin(); i != tags.end();
                ++i ) {
                aux_outdegree[ *i ] += 1;
            }
            tags.clear();
        } else if( Chromosome::IsInteraction()( *i ) ) {
            int ia = dynamic_cast< Interaction* >( *i )->tag();
            // bsite cp number
            aux_bsites[ ia ] += 1;
            // unique interactions per gene
            tags.insert( ia );
        }
#ifdef DOUBLESTRANDBREAKS
          else if( Chromosome::IsRepeat()( *i ) ) {
            ++nr_repeats_;
            /*tags.clear();*/
        } else if( Chromosome::IsRetroposon()( *i ) ) {
            ++nr_retroposons_;
        }
#endif
    }
    // and add final histos
    addHisto( nr_genes_, aux_genes );
    addHisto( nr_bsites_, aux_bsites );
    addHisto( in_degree_, aux_indegree );
    addHisto( out_degree_, aux_outdegree );
}

void
dorsalfin::LogGenomeStats::addHisto( genomehisto &gh, const histo &h ) {
    for( histo::const_iterator i( h.begin() ); i != h.end(); ++i ) {
        ++gh[ i->first ][ i->second ];
    }
}

//
// Child observer of GenomeStats
//
dorsalfin::LogCsvGenes::LogCsvGenes( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# copy nr range, individuals per gene\n";
}

void
dorsalfin::LogCsvGenes::update( Subject *s ) {
    // write histogram of each gene cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    LogGenomeStats::genomehisto gh = lgs->nrGenes();
    // determine min/max
    int min = lgs->sample(), max = 0;
    for( LogGenomeStats::genomehisto::const_iterator i( gh.begin() );
        i != gh.end(); ++i ) {
        if( min > i->second.begin()->first ) min = i->second.begin()->first;
        if( max < i->second.rbegin()->first ) max = i->second.rbegin()->first;
    }
    // output (assuming 
    int nn = Network::referenceIds().size();
    for( int c = min; c != max+1; ++c ) {
        *log_ << c << "\t";
        for( int i = 0; i != nn; ++i ) {
            // does this work?
            *log_ << gh[ i ][ c ] << "\t";
        }
        *log_ << "\n";
    }
    *log_ << "\n\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 2nd child observer of GenomeStats
//
dorsalfin::LogCsvBindingSites::LogCsvBindingSites( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# cp nr range, individuals per binding site\n";
}

void
dorsalfin::LogCsvBindingSites::update( Subject *s ) {
    // write histogram of each gene cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    LogGenomeStats::genomehisto gh = lgs->nrBindingSites();
    // determine min/max
    int min = lgs->sample(), max = 0;
    for( LogGenomeStats::genomehisto::const_iterator i( gh.begin() );
        i != gh.end(); ++i ) {
        if( min > i->second.begin()->first ) min = i->second.begin()->first;
        if( max < i->second.rbegin()->first ) max = i->second.rbegin()->first;
    }
    // output (assuming 
    int nn = Network::referenceIds().size();
    for( int c = min; c != max+1; ++c ) {
        *log_ << c << "\t";
        for( int i = 0; i != nn; ++i ) {
            // does this work?
            *log_ << gh[ i ][ c ] << "\t";
       }
        *log_ << "\n";
    }
    *log_ << "\n\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 4th child observer of GenomeStats
//
dorsalfin::LogCsvInDegree::LogCsvInDegree( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# range, individuals per nr of unique bsites per gene\n";
}

void
dorsalfin::LogCsvInDegree::update( Subject *s ) {
    // write histogram of each gene cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    LogGenomeStats::genomehisto gh = lgs->inDegree();
    // determine min/max
    int min = lgs->sample(), max = 0;
    for( LogGenomeStats::genomehisto::const_iterator i( gh.begin() );
        i != gh.end(); ++i ) {
        if( min > i->second.begin()->first ) min = i->second.begin()->first;
        if( max < i->second.rbegin()->first ) max = i->second.rbegin()->first;
    }
    // output (assuming 
    int nn = Network::referenceIds().size();
    for( int c = min; c != max+1; ++c ) {
        *log_ << c << "\t";
        for( int i = 0; i != nn; ++i ) {
            // does this work?
            *log_ << gh[ i ][ c ] << "\t";
        }
        *log_ << "\n";
    }
    *log_ << "\n\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 5th child observer of GenomeStats
//
dorsalfin::LogCsvOutDegree::LogCsvOutDegree( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# range, individuals per nr of unique targets per gene\n";
}

void
dorsalfin::LogCsvOutDegree::update( Subject *s ) {
    // write histogram of each gene cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    LogGenomeStats::genomehisto gh = lgs->outDegree();
    // determine min/max
    int min = lgs->sample(), max = 0;
    for( LogGenomeStats::genomehisto::const_iterator i( gh.begin() );
        i != gh.end(); ++i ) {
        if( min > i->second.begin()->first ) min = i->second.begin()->first;
        if( max < i->second.rbegin()->first ) max = i->second.rbegin()->first;
    }
    // output (assuming 
    int nn = Network::referenceIds().size();
    for( int c = min; c != max+1; ++c ) {
        *log_ << c << "\t";
        for( int i = 0; i != nn; ++i ) {
            // does this work?
            *log_ << gh[ i ][ c ] << "\t";
        }
        *log_ << "\n";
    }
    *log_ << "\n\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 7th child observer of GenomeStats
//
dorsalfin::LogCsvRetroposons::LogCsvRetroposons( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# retroposons repeats\n";
}

void
dorsalfin::LogCsvRetroposons::update( Subject *s ) {
    // write averages of each gene's reference state 
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    *log_ << lgs->nrRetroposons() / cc << "\t" << lgs->nrRepeats() / cc << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

