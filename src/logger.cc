//
// Implementation of a few observers that dump their data in files.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "logger.hh"

#include "environment.hh"
#include "population.hh"
#include "genreg_agent.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "retroposon.hh"
#include "gene.hh"

#include "duo_agent.hh"

//
// Counting double strand breaks and other mutations
//
dorsalfin::LogCsvMutations::LogCsvMutations( std::string fname ) 
    : AsyncLogObserver( /*s*/ ), time_( -1 ),
     cps_( 4, 0 ), rms_( 4, 0 ), thrs_( 4, 0 ), iacps_( 4, 0 ), 
     iarms_( 4, 0 ), ians_( 4, 0 ), iats_( 4, 0 ), iaws_( 4, 0 ) {
    openLog( fname );
    *log_ << "# WARNING Only handles one module\n";
    *log_ << "# col 1: time\n";
    *log_ << "# col 2-5: pos, neg, neu (cp)\n";
    *log_ << "# col 6-9: pos, neg, neu (rm)\n";
    *log_ << "# col 10-13: pos, neg, neu (thr)\n";
    *log_ << "# col 14-17: pos, neg, neu (ia cp)\n";
    *log_ << "# col 18-21: pos, neg, neu (ia rm)\n";
    *log_ << "# col 22-25: pos, neg, neu (ia nw)\n";
    *log_ << "# col 26-29: pos, neg, neu (ia tag)\n";
    *log_ << "# col 30-34: pos, neg, neu (ia weight)\n";
#ifdef DOUBLESTRANDBREAKS
    *log_ << "# col 35-39: pos, neg, neu (dsb) TODO\n";
#endif
}

void
dorsalfin::LogCsvMutations::write() {
    if( std::accumulate( cps_.begin(), cps_.end(), 0u ) > 0 or
        std::accumulate( rms_.begin(), rms_.end(), 0u ) > 0 or
        std::accumulate( thrs_.begin(), thrs_.end(), 0u ) > 0 or 
        std::accumulate( iacps_.begin(), iacps_.end(), 0u ) > 0 or
        std::accumulate( iarms_.begin(), iarms_.end(), 0u ) > 0 or
        std::accumulate( ians_.begin(), ians_.end(), 0u ) > 0 or
        std::accumulate( iats_.begin(), iats_.end(), 0u ) > 0 or
        std::accumulate( iaws_.begin(), iaws_.end(), 0u ) > 0 ) {
        // new timestep, write old stuff if nonzero
        *log_ << time_ << "\t\t";
        std::copy( cps_.begin(), cps_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( rms_.begin(), rms_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( thrs_.begin(), thrs_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iacps_.begin(), iacps_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iarms_.begin(), iarms_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( ians_.begin(), ians_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iats_.begin(), iats_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iaws_.begin(), iaws_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\n";
        
        // do not want accumulative: reset to zero
        std::fill( cps_.begin(), cps_.end(), 0 );
        std::fill( rms_.begin(), rms_.end(), 0 );
        std::fill( thrs_.begin(), thrs_.end(), 0 );
        std::fill( iacps_.begin(), iacps_.end(), 0 );
        std::fill( iarms_.begin(), iarms_.end(), 0 );
        std::fill( ians_.begin(), ians_.end(), 0 );
        std::fill( iats_.begin(), iats_.end(), 0 );
        std::fill( iaws_.begin(), iaws_.end(), 0 );
    }
#ifdef DEBUG
    log_->flush();
#endif
}    

void
dorsalfin::LogCsvMutations::update( Subject *s ) {
    // We are receiving two agents..
    DuoAgent *da = dynamic_cast< DuoAgent* >( s );
    GenRegAgent *m1 = dynamic_cast< GenRegAgent* >( da->first );
    GenRegAgent *m2 = dynamic_cast< GenRegAgent* >( da->second );
    if( m1 and m2 ) {
        // get the time
        if( m1->myTag().time > time_ ) {
            write();
            time_ = m1->myTag().time;
        }
        
        understand( *m1 );
        understand( *m2 );
   }
}

void
dorsalfin::LogCsvMutations::understand( const GenRegAgent &n ) {
    std::vector< uint > mm( n.nrMutations() );
    // gene dist differences
    double dd = n.deltaDistance();
    int ds = n.deltaSize();

    if( dd < 0 ) {
        // assume something positive happened
        cps_[ POS ] += mm[ Chromosome::CP_G ];
        rms_[ POS ] += mm[ Chromosome::RM_G ];
        thrs_[ POS ] += mm[ Chromosome::THR ];
        iacps_[ POS ] += mm[ Chromosome::CP_IA ];
        iarms_[ POS ] += mm[ Chromosome::RM_IA ];
        ians_[ POS ] += mm[ Chromosome::NW_IA ];
        iats_[ POS ] += mm[ Chromosome::T_IA ];
        iaws_[ POS ] += mm[ Chromosome::W_IA ];
    } else if( dd > 0 ) {
        // negative
        cps_[ NEG ] += mm[ Chromosome::CP_G ];
        rms_[ NEG ] += mm[ Chromosome::RM_G ];
        thrs_[ NEG ] += mm[ Chromosome::THR ];
        iacps_[ NEG ] += mm[ Chromosome::CP_IA ];
        iarms_[ NEG ] += mm[ Chromosome::RM_IA ];
        ians_[ NEG ] += mm[ Chromosome::NW_IA ];
        iats_[ NEG ] += mm[ Chromosome::T_IA ];
        iaws_[ NEG ] += mm[ Chromosome::W_IA ];
    } else {
        if( ds > 0 ) {
            // neutral
            cps_[ NEU ] += mm[ Chromosome::CP_G ];
            iacps_[ NEU ] += mm[ Chromosome::CP_IA ];
            ians_[ NEU ] += mm[ Chromosome::NW_IA ];
        } else if( ds < 0 ) {
            rms_[ NEU ] += mm[ Chromosome::RM_G ];
            iarms_[ NEU ] += mm[ Chromosome::RM_IA ];
        } else {
            thrs_[ NEU ] += mm[ Chromosome::THR ];
            iats_[ NEU ] += mm[ Chromosome::T_IA ];
            iaws_[ NEU ] += mm[ Chromosome::W_IA ];
        }
    }
}

//
// Simple csv observer for the distances
//
dorsalfin::LogCsvDistances::LogCsvDistances( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvDistances::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the distances
    std::vector< double > distances;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        if( not na->justBorn() ) {
            distances.push_back( na->distance() );
        }
    }
    
    // calculate min, median, mean and variance
    if( not distances.empty() ) {
        std::stable_sort( distances.begin(), distances.end() ); 
        double max = distances.back();
        double avg = mean( distances.begin(), distances.end() ); 
        
        // median
        double median = 0.0;
        uint n = distances.size();
        if( n % 2 == 0 ) {
            median = ( distances[ n / 2 ] + distances[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = distances[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( distances.begin(), distances.end(), avg ) );
        
        *log_ << max << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
    } else {
        *log_ << "# Empty grid...\n";
    }
//#ifdef DEBUG
    log_->flush();
//#endif
}

void
dorsalfin::LogCsvDistances::writeHeader() {
    *log_ << "# max, median, mean, stdev of genotypic distances\n";
}


//
// Simple csv observer for the distances in a histogram
//
dorsalfin::LogCsvHistoDistances::LogCsvHistoDistances( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvHistoDistances::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the distances
    std::map< uint, uint > distances;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        if( not na->justBorn() ) {
            ++distances[ na->distance() ];
        }
    }
    
    // output
    ulong time = pop->model()->now();
    for( std::map< uint, uint >::iterator i( distances.begin() );
        i != distances.end(); ++i ) {
        *log_ << time << "\t" << i->first << "\t" << i->second << "\n";
    }
    
//#ifdef DEBUG
    log_->flush();
//#endif
}

void
dorsalfin::LogCsvHistoDistances::writeHeader() {
    *log_ << "# time, genotypic distance, nr of agents\n";
}


//
// Simple csv observer for the lengths in a histogram
//
dorsalfin::LogCsvHistoLengths::LogCsvHistoLengths( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvHistoLengths::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the lengths
    std::map< uint, uint > lens;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        //if( !na->justBorn() ) {
            ++lens[ na->length() ];
        //}
    }
    
    // output
    ulong time = pop->model()->now();
    for( std::map< uint, uint >::iterator i( lens.begin() );
        i != lens.end(); ++i ) {
        *log_ << time << "\t" << i->first << "\t" << i->second << "\n";
    }
    
//#ifdef DEBUG
    log_->flush();
//#endif
}

void
dorsalfin::LogCsvHistoLengths::writeHeader() {
    *log_ << "# time, genome length, nr of agents\n";
}


//
// Simple csv observer for the scores
//
dorsalfin::LogCsvScores::LogCsvScores( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvScores::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the scores
    std::vector< double > scores;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        scores.push_back( i->first->score() );
    }
    
    // calculate max, median, mean and variance
    if( !scores.empty() ) {
        std::stable_sort( scores.begin(), scores.end() ); 
        double max = scores.back();
        double avg = mean( scores.begin(), scores.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores.size();
        if( n % 2 == 0 ) {
            median = ( scores[ n / 2 ] + scores[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( scores.begin(), scores.end(), avg ) );
       
        *log_ << max << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
    } else {
        *log_ << "# Empty grid...\n";
    }
    log_->flush();
}

void
dorsalfin::LogCsvScores::writeHeader() {
    *log_ << "# max, median, mean, variance of raw fitness scores\n";
}

//
// Another class
// 
dorsalfin::LogXmlGenomes::LogXmlGenomes( std::string dname, long i ) 
    : LogObserver( i ) {
    // create dir with given name
    // within dir, use some logical name
    // f.i. timestep.dot
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
dorsalfin::LogXmlGenomes::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    *log_ << "<generation time=\"" << pop->generation() << "\">\n"
          << *pop << "</generation>\n";
    closeLog();
}

void
dorsalfin::LogXmlGenomes::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation dorsal_version=\"" << VERSION << "\">\n";
}

void
dorsalfin::LogXmlGenomes::writeFooter() {
    *log_ << "</simulation>\n";
}

std::string
dorsalfin::LogXmlGenomes::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".xml";
    return result.str();
}

//
// Another class
// 
dorsalfin::LogXmlSampleGenomes::LogXmlSampleGenomes( std::string dname, long i ) 
    : LogObserver( i ), fraction_( 0.1 ), sample_generator_( 2 ), 
      sample_uniform_( sample_generator_ ) {
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
dorsalfin::LogXmlSampleGenomes::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();
   
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    *log_ << "<generation time=\"" << pop->generation() << "\">\n";
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        if( sample_uniform_() < fraction_ ) {
            *log_ << *( i->first );
        }
    }
    *log_ << "</generation>\n";
    closeLog();
}

void
dorsalfin::LogXmlSampleGenomes::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation dorsal_version=\"" << VERSION << "\">\n";
}

void
dorsalfin::LogXmlSampleGenomes::writeFooter() {
    *log_ << "</simulation>\n";
}

std::string
dorsalfin::LogXmlSampleGenomes::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".xml";
    return result.str();
}


//
// Log the size of the population
//
dorsalfin::LogPopulationSize::LogPopulationSize(
    std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogPopulationSize::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    *log_ << pop->nrAgents() << std::endl;
}

void
dorsalfin::LogPopulationSize::writeHeader() {
    *log_ << "# population size\n";
}


//
// Log the output of the agents
//
dorsalfin::LogCsvPopulationBitset::LogCsvPopulationBitset(
    std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvPopulationBitset::initialize( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    BitEnvironment *be = dynamic_cast< BitEnvironment * >( 
        &pop->model()->environment() );
    // get bits and rotate them all
    std::vector< boost::dynamic_bitset<> > bits = be->bits();
    for( std::vector< boost::dynamic_bitset<> >::iterator j( bits.begin() );
        j != bits.end(); ++j ) {
        uint offset = std::distance( bits.begin(), j );
        for( uint i = 0; i != j->size(); ++i ) {
            rotated_[ *j ] = offset + i;
            rotate( *j, 1 );
        }
    }
}

void
dorsalfin::LogCsvPopulationBitset::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();
    
    // which environment?
    NoRotateBitEnvironment *be = dynamic_cast< NoRotateBitEnvironment * >( 
        &pop->model()->environment() );
    
    std::map< std::pair< uint, uint >, uint > outputs;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const EnvAgent *na = dynamic_cast< EnvAgent* >( i->first );
        if( na->reproduced() ) {
            if( be ) {
                ++outputs[ std::make_pair( 0, na->output().size() ) ];
            } else {
                ++outputs[ std::make_pair( rotated_[ na->input() ],
                    na->output().size() ) ];
            }                
        }
    }

    ulong time = pop->model()->now();
    std::map< std::pair< uint, uint >, uint >::iterator i;
    for( i = outputs.begin(); i != outputs.end(); ++i ) {
        *log_ << time << "\t" 
            << i->first.first << "\t" << i->first.second << "\t" 
            << i->second << "\n";
    }
    *log_ << "\n\n";
    
//#ifdef DEBUG
    log_->flush();
//#endif
}

void
dorsalfin::LogCsvPopulationBitset::writeHeader() {
    *log_ << "# histo of output sequences: start, length, count\n";
}


//
// Log what's in the environment and how often
//
dorsalfin::LogCsvEnvironmentBitset::LogCsvEnvironmentBitset(
    std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvEnvironmentBitset::doUpdate( Subject *s ) {
    BitEnvironment *env = dynamic_cast< BitEnvironment * >( s );
    BitEnvironment::bitset_grid grid = env->grid();

    // get maximum length of a bitstring
    uint maxlen = env->bits().front().size();
    // and go!
    uint empty = 0;    
    std::map< std::string, uint > outputs;
    int n = grid.shape()[ 0 ];
    int m = grid.shape()[ 1 ];
    for( int i = 0; i != n; ++i ) {
        for( int j = 0; j != m; ++j ) {
            std::string aux;
            boost::to_string( grid[ i ][ j ], aux );
            aux.insert( aux.begin(), maxlen - aux.size(), '.' );
            outputs[ aux ] += 1;
        }
    }

    *log_ << "# empty\t" << empty << "\n";
    ulong time = env->model()->now();
    std::map< std::string, uint >::iterator i;
    for( i = outputs.begin(); i != outputs.end(); ++i ) {
        *log_ << time << "\t" << i->first << "\t" << i->second << "\n";
    }
    *log_ << "\n\n";
    
//#ifdef DEBUG
    log_->flush();
//#endif
}

void
dorsalfin::LogCsvEnvironmentBitset::writeHeader() {
    *log_ << "# histo of environmental sequences\n";
}

