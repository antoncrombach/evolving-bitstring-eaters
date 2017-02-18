//
// Some observers of the population.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_LOGGENOME_H_
#define _DORSALFIN_LOGGENOME_H_

#include "defs.hh"
#include "observer.hh"
#include "population.hh"
#include "stream_manager.hh"

namespace dorsalfin {

    /// \class LogGenomeStats
    ///
    /// Gather a bunch of stats and counts from each genome and write to 
    /// some different files different aspects of these genomes.
    ///
    /// How to do this elegantly? Multiple inheritance, this observer
    /// is being observed by three others. Those three are the ones actually
    /// writing to file. This observer notifies them if they need to write...
    class LogGenomeStats : public LogObserver, public Subject {
        public:
        typedef std::map< int, int > histo;
        typedef std::map< int, histo > genomehisto;
        
        public:
        LogGenomeStats( long );
        virtual ~LogGenomeStats() {}
        
        virtual void initialize( Subject * ) {}
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        const genomehisto & nrGenes() const;
        const genomehisto & nrBindingSites() const;
        const genomehisto & inDegree() const;
        const genomehisto & outDegree() const;
        
        int nrRepeats() const;
        int nrRetroposons() const;
        
        int sample() const;
        
        private:
        void doCounting( Chromosome * );
        void addHisto( genomehisto &, const histo & );
        
        private:
        int sample_, nr_repeats_, nr_retroposons_;
        genomehisto nr_genes_, nr_bsites_, in_degree_, out_degree_;
    };

    inline const LogGenomeStats::genomehisto &
    LogGenomeStats::nrGenes() const 
    { return nr_genes_; }

    inline const LogGenomeStats::genomehisto &
    LogGenomeStats::nrBindingSites() const 
    { return nr_bsites_; }

    inline const LogGenomeStats::genomehisto &
    LogGenomeStats::inDegree() const 
    { return in_degree_; }

    inline const LogGenomeStats::genomehisto & 
    LogGenomeStats::outDegree() const
    { return out_degree_; }

    inline int LogGenomeStats::nrRepeats() const
    { return nr_repeats_; }

    inline int LogGenomeStats::nrRetroposons() const 
    { return nr_retroposons_; }

    inline int LogGenomeStats::sample() const 
    { return sample_; }    

    
    class LogCsvGenes : public AsyncLogObserver {
        public:
        LogCsvGenes( std::string );
        virtual ~LogCsvGenes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvBindingSites : public AsyncLogObserver {
        public:
        LogCsvBindingSites( std::string );
        virtual ~LogCsvBindingSites() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };
    
    class LogCsvInDegree : public AsyncLogObserver {
        public:
        LogCsvInDegree( std::string );
        virtual ~LogCsvInDegree() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvOutDegree : public AsyncLogObserver {
        public:
        LogCsvOutDegree( std::string );
        virtual ~LogCsvOutDegree() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvRetroposons : public AsyncLogObserver {
        public:
        LogCsvRetroposons( std::string );
        virtual ~LogCsvRetroposons() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };
}
#endif

