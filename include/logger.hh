//
// Some observers of the population.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_LOGGER_H_
#define _DORSALFIN_LOGGER_H_

#include "defs.hh"
#include "observer.hh"
#include "population.hh"
#include "stream_manager.hh"

namespace dorsalfin {

    /// \class LogCsvMutations
    ///
    /// Log the mutations of a timestep. If no mutations, nothing is logged.
    class LogCsvMutations : public AsyncLogObserver {
        public:
        enum clsf { POS, NEG, NEU, UNK };
        
        public:
        LogCsvMutations( std::string );
        virtual ~LogCsvMutations() {}

        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void write();
        void understand( const GenRegAgent & );
        
        private:
        long time_;
        std::vector< uint > cps_, rms_, thrs_, 
            iacps_, iarms_, ians_, iats_, iaws_;
    };    

    /// \class LogCsvScores
    ///
    /// Log min, mean, median and stddev of fitness scores.
    class LogCsvScores : public LogObserver {
        public:
        LogCsvScores( std::string, long );
        virtual ~LogCsvScores() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvDistances
    ///
    /// Log mean, min, median and stddev of distances.
    class LogCsvDistances : public LogObserver {
        public:
        LogCsvDistances( std::string, long );
        virtual ~LogCsvDistances() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvHistoDistances
    ///
    /// Log histogram of distances.
    class LogCsvHistoDistances : public LogObserver {
        public:
        LogCsvHistoDistances( std::string, long );
        virtual ~LogCsvHistoDistances() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvHistoLengths
    ///
    /// Log histogram of genome lengths
    class LogCsvHistoLengths : public LogObserver {
        public:
        LogCsvHistoLengths( std::string, long );
        virtual ~LogCsvHistoLengths() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };    
    
    /// \class LogCsvEnvironmentBitset
    ///
    /// Log the input of individuals
    class LogCsvEnvironmentBitset : public LogObserver {
        public:
        LogCsvEnvironmentBitset( std::string, long );
        virtual ~LogCsvEnvironmentBitset() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvPopulationBitset
    ///
    /// Log the output of individuals
    class LogCsvPopulationBitset : public LogObserver {
        public:
        LogCsvPopulationBitset( std::string, long );
        virtual ~LogCsvPopulationBitset() {}
        
        virtual void initialize( Subject * );
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        
        private:
        std::map< boost::dynamic_bitset<>, uint > rotated_;
    };
    
    /// \class LogXmlGenomes
    ///
    /// Log entire populations.. at least the genomes.
    class LogXmlGenomes : public LogObserver {
        public:
        LogXmlGenomes( std::string, long );
        virtual ~LogXmlGenomes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize()
        { writeFooter(); }
        
        private:
        void writeHeader();
        void writeFooter();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

    /// \class LogXmlSampleGenomes
    ///
    /// Log a few individuals sampled from populations.
    class LogXmlSampleGenomes : public LogObserver {
        public:
        LogXmlSampleGenomes( std::string, long );
        virtual ~LogXmlSampleGenomes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize()
        { writeFooter(); }
        
        private:
        void writeHeader();
        void writeFooter();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
        
        double fraction_;
        base_generator_type sample_generator_;
        uniform_gen_type sample_uniform_;
    };

    /// \class LogPopulationSize
    ///
    /// Log the number of individuals.
    class LogPopulationSize : public LogObserver {
        public:
        LogPopulationSize( std::string, long );
        virtual ~LogPopulationSize() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };
}
#endif

