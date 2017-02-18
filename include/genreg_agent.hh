//
// Interface of a simple agent with a genome and network.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_GENREGAGENT_H_
#define _DORSALFIN_GENREGAGENT_H_

#include "defs.hh"
#include "agent.hh"
#include "genome.hh"
#include "network.hh"

namespace dorsalfin {

    /// \class GenRegAgent
    /// \brief Slowly including genome and network.
    class GenRegAgent : public Agent {
        public:
        struct ParentInfo {
            int size_;
            int distance_;
        };
        
        public:
        /// Init
        virtual void initialise();

        /// Set the genome (rebuild network?)
        void genome( Genome * );
        /// Get the genome
        const Genome* genome() const;
        /// Get the network
        const Network* network() const;

        /// Get the score of a \c GenRegAgent
        virtual double score() const;
        /// Get genotypical distance to target
        int distance() const;
        /// Get genome length
        int length() const;
        int size() const;

        /// Set if we read the thing from file
        void readFromFile( bool );
        
        /// Get a few mutation counters
        std::vector< uint > nrMutations() const;
        /// Get the nr of state changes
        uint nrPerturbations() const;
        /// Get size difference with parent
        int deltaSize() const;
        /// Get genotype distance difference with parent
        int deltaDistance() const;

        protected:
        /// Default ctor.
        GenRegAgent();
        /// Ctor.
        GenRegAgent( uint );
        /// Ctor.
        GenRegAgent( uint, Genome *g );
        /// Copy ctor.
        GenRegAgent( const GenRegAgent & );
        /// Destructor.
        virtual ~GenRegAgent();
        /// Cloning an agent.
        /*virtual Agent* clone() const;*/
        /// Copy the agent into \c this
        virtual void copy( const Agent & );
        
        public:
        /// Set death rate of agent
        static void deathRate( double );
        /// Get death rate of agent
        static double deathRate();
        /// Set maximum distance
        static void maxDistance( int );
        /// Get maximum distance
        static int maxDistance();
        /// Set max genome size, after which a penalty is applied
        static void maxGenomeSize( int );
        /// Set max # retroposons, after which a penalty is applied
        static void maxTposons( int );
        /// Set penalty per unit over the max size/number
        static void genomePenaltyRate( double );
        /// Set penalty per retroposon too much
        static void retroposonPenaltyRate( double );
        /// Set initial state of network
        static void initialState( const boost::dynamic_bitset<> & );
            
        protected:
        /// Raw score aka distance to attractor
        int distance_;
        /// Genome
        Genome *genome_;
        /// Regulatory network
        Network *network_;
        /// Parent info
        ParentInfo parent_;
        /// Flag for network init
        bool read_from_file_;
        
        protected:
        /// Death rate
        static double death_rate_;
        /// Maximum distance we look at
        static int max_dist_;
        /// Penalty if genome size is greater
        static int penalty_genome_;
        /// Penalty if more retroposons than allowed
        static int penalty_tposons_;
        /// Penalty coefficient for genome
        static double penalty_genome_rate_;
        /// Penalty coefficient for retroposon
        static double penalty_tposons_rate_;
        /// Initial state of network
        static boost::dynamic_bitset<> init_state_;
    };

    inline const Genome* GenRegAgent::genome() const
    { return genome_; }
    
    inline const Network* GenRegAgent::network() const
    { return network_; }
    
    inline int GenRegAgent::distance() const
    { return distance_; }

    inline int GenRegAgent::length() const
    { return genome_->length(); }
    
    inline int GenRegAgent::size() const
    { return genome_->length(); }

    inline int GenRegAgent::deltaDistance() const
    { return distance_ - parent_.distance_; }
    
    inline int GenRegAgent::deltaSize() const
    { return genome_->length() - parent_.size_; }

    inline std::vector< uint > GenRegAgent::nrMutations() const
    { return genome_->nrMutations(); }
    
    inline uint GenRegAgent::nrPerturbations() const
    { return genome_->nrPerturbations(); }
    
    inline void GenRegAgent::readFromFile( bool b )
    { read_from_file_ = b; }
}
#endif

