//
// Implementation of a few observers that dump their data in files.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "loggrid.hh"

#include "environment.hh"
#include "population.hh"
#include "genreg_agent.hh"
//#include "genome.hh"
//#include "chromosome.hh"
//#include "retroposon.hh"
//#include "gene.hh"

//#include "duo_agent.hh"

//
// Simple csv observer for the environment
//
dorsalfin::LogCsvEnvironmentGrid::LogCsvEnvironmentGrid( 
        std::string dname, long i ) : LogObserver( i ) {
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
dorsalfin::LogCsvEnvironmentGrid::doUpdate( Subject *s ) {
    BitEnvironment *env = dynamic_cast< BitEnvironment * >( s );
    BitEnvironment::bitset_grid grid = env->grid();
    // begin of file
    openLog( unique_name( env->model()->now() ) );
    writeHeader();
    // and a label of the class (not portable between compilers!)
    std::string outname( typeid( *s ).name() );
    outname.erase( 0, outname.find( "19" ) + 2 );
    outname.erase( outname.size() - 1, 1 );
    *log_ << outname << "\t";
    // with the dimensions of the matrix to follow
    int n = grid.shape()[ 0 ];
    int m = grid.shape()[ 1 ];
    *log_ << n << "\t" << m << "\n";
    for( int i = 0; i != n; ++i ) {
        for( int j = 0; j != m; ++j ) {
            *log_ << "0b" << grid[ i ][ j ] << "\t";        
        }
        *log_ << "\n";
    }
    closeLog();
}

void
dorsalfin::LogCsvEnvironmentGrid::writeHeader() {
    *log_ << "# time, matrix of environment bits\n";
}

std::string
dorsalfin::LogCsvEnvironmentGrid::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".csv";
    return result.str();
}

//
// And now a population dumper, at least it dumps a few specific feats
//
dorsalfin::LogCsvPopulationGrid::LogCsvPopulationGrid( std::string dname, long i ) 
    : LogObserver( i ) {
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
dorsalfin::LogCsvPopulationGrid::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_grid grid = pop->grid();
    
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    // and a label of the feat
    *log_ << "gene distance" << "\t";
    // with the dimensions of the matrix to follow
    int n = grid.shape()[ 0 ];
    int m = grid.shape()[ 1 ];
    *log_ << n << "\t" << m << "\n";
    for( int i = 0; i != n; ++i ) {
        for( int j = 0; j != m; ++j ) {
            if( grid[ i ][ j ] != 0 ) {
                if( grid[ i ][ j ]->justBorn() ) {
                    *log_ << 0.0 << "\t";
                } else {
                    *log_ << grid[ i ][ j ]->score() << "\t";
                }
            } else {
                *log_ << -1 << "\t";
            }
        }
        *log_ << "\n";
    }
    
    // we want to see the output too
    *log_ << "\n\noutput bits" << "\t";
    // with the dimensions of the matrix to follow
    *log_ << n << "\t" << m << "\n";
    for( int i = 0; i != n; ++i ) {
        for( int j = 0; j != m; ++j ) {
            if( grid[ i ][ j ] != 0 ) {
                const EnvAgent *na = 
                    dynamic_cast< EnvAgent* >( grid[ i ][ j ] );
                if( na->justBorn() ) {
                    *log_ << -2 << "\t";
                } else {
                    *log_ << "0b" << na->output() << "\t";
                }
            } else {
                *log_ << -1 << "\t";
            }
        }
        *log_ << "\n";
    }
    
    closeLog();
}

void
dorsalfin::LogCsvPopulationGrid::writeHeader() {
    *log_ << "# grid size and features in matrix format,";
    *log_ << " separate by newlines\n";
}

std::string
dorsalfin::LogCsvPopulationGrid::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".csv";
    return result.str();
}


//
// Record only a single line in the middle of the population grid
//
dorsalfin::LogCsvPopGridLine::LogCsvPopGridLine( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvPopGridLine::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_grid grid = pop->grid();
    
    // the line and its length
    int n = grid.shape()[ 0 ] / 2;
    int m = grid.shape()[ 1 ];
    for( int j = 0; j != m; ++j ) {
        if( grid[ n ][ j ] != 0 ) {
            const EnvAgent *na = 
                dynamic_cast< EnvAgent* >( grid[ n ][ j ] );
            if( na->justBorn() ) {
                *log_ << -2 << "\t";
            } else {
                *log_ << "0b" << na->output() << "\t";
            }
        } else {
            *log_ << -1 << "\t";
        }
    }
    *log_ << "\n";
}

void
dorsalfin::LogCsvPopGridLine::writeHeader() {
    *log_ << "# output of agents on the middle line\n";
}

//
// Record only a single line in the middle of the environment grid
//
dorsalfin::LogCsvEnvGridLine::LogCsvEnvGridLine( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
dorsalfin::LogCsvEnvGridLine::doUpdate( Subject *s ) {
    BitEnvironment *env = dynamic_cast< BitEnvironment * >( s );
    BitEnvironment::bitset_grid grid = env->grid();
    
    // the line and its length
    int n = grid.shape()[ 0 ] / 2;
    int m = grid.shape()[ 1 ];
    for( int j = 0; j != m; ++j ) {
        *log_ << "0b" << grid[ n ][ j ] << "\t";
    }
    *log_ << "\n";
}

void
dorsalfin::LogCsvEnvGridLine::writeHeader() {
    *log_ << "# environment on the middle line\n";
}

