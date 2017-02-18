//
// Some observers of the population.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_LOGGRID_H_
#define _DORSALFIN_LOGGRID_H_

#include "defs.hh"
#include "observer.hh"
#include "population.hh"
#include "stream_manager.hh"

namespace dorsalfin {

    /// \class LogCsvEnvironmentGrid
    ///
    /// Log change of environment.
    class LogCsvEnvironmentGrid : public LogObserver {
        public:
        LogCsvEnvironmentGrid( std::string, long );
        virtual ~LogCsvEnvironmentGrid() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}

        private:
        void writeHeader();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };        

    /// \class LogCsvPopulationGrid
    ///
    /// Log the fitness in a matrix format. The resulting file can be used
    /// to generate movies of the grid.
    class LogCsvPopulationGrid : public LogObserver {
        public:
        LogCsvPopulationGrid( std::string, long );
        virtual ~LogCsvPopulationGrid() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

    /// \class LogCsvEnvGridLine
    ///
    /// Log bitsets of environment at middle grid line.
    class LogCsvEnvGridLine : public LogObserver {
        public:
        LogCsvEnvGridLine( std::string, long );
        virtual ~LogCsvEnvGridLine() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}

        private:
        void writeHeader();
    };        

    /// \class LogCsvPopGridLine
    ///
    /// Log output of agents at middle grid line.
    class LogCsvPopGridLine : public LogObserver {
        public:
        LogCsvPopGridLine( std::string, long );
        virtual ~LogCsvPopGridLine() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}

        private:
        void writeHeader();
    };
}
#endif

