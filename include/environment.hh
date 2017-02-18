//
// Environment base class and childs.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_ENVIRONMENT_H_
#define _DORSALFIN_ENVIRONMENT_H_

#include "defs.hh"
#include "subject.hh"
#include "population.hh"
#include "agent.hh"
#include "env_agent.hh"

namespace dorsalfin {

    /// \class Environment
    /// \brief Abstract base class for modelling the environment.
    class Environment : public Subject {
        public:
        /// Destructor (empty).
        virtual ~Environment() {};
        
        /// Initialise
        virtual void initialise( 
            const std::vector< boost::dynamic_bitset<> > & ) = 0;
        /// The environment changes, the population is alerted
        virtual void fluctuate( long ) = 0;
        /// The environment is changed by something outside
        virtual void update( const Agent *, const location & ) = 0;
        /// Finish of the environment
        virtual void finish();
        
        /// Set a reference to the model
        void model( Model * );
        /// Get a reference to the model
        Model* model() const;
        
        protected:
        /// Hidden constructor
        Environment( uint );
        
        protected:
        /// The environment has its own uniform random nr generator
        base_generator_type generator_env_;
        uniform_gen_type uniform_env_;
        
        /// Reference to the model
        Model *model_;
    };
    
    inline void Environment::model( Model *m )
    { model_ = m; }
    
    inline Model* Environment::model() const
    { return model_; }
    
    /// \class BitEnvironment
    /// \brief Environment with bitstrings
    ///
    /// Bitstrings are changed by individuals...
    class BitEnvironment : public Environment {
        public:
        /// Typedef of the bitstring grid
        typedef boost::multi_array< boost::dynamic_bitset<>, 2 > bitset_grid;
        /// Range of grid
        typedef bitset_grid::index_range range;

        public:
        /// Destructor
        ~BitEnvironment() {};
        
        /// Initialise grid with bitset
        virtual void initialise( 
            const std::vector< boost::dynamic_bitset<> > & );
        
        /// Get the bitset from a certain location
        void getBits( boost::dynamic_bitset<> &, const location & ) const;
        /// Set the bitset at a certain location
        void setBits( const boost::dynamic_bitset<> &, const location & );
        
        /// Get the moore neighbourhood of the given location (from shadow)
        void moore( std::vector< boost::dynamic_bitset<> > &,
            const location & ) const;
        /// Get grid (observer needs it)
        const bitset_grid & grid() const;
        /// Get original bitsets
        const std::vector< boost::dynamic_bitset<> > & bits() const;

        /// Margolus diffusion
        void diffusion();
        
        /// Read environmental grid
        void readGrid( std::ifstream & );
        
        public:
        /// Set diffusion period
        static void setDiffusion( double );
        
        protected:
        /// Hidden constructor
        BitEnvironment() : Environment( 18 ), rw_grid_(), focus_() {};
        /// Constructor with dimensions
        BitEnvironment( int, int, uint );

        private:
        // Auxilary function for Margolis diffusion
        void rotate( uint, uint, uint, uint );
        
        protected:
        // As long as we read/write only to empty cells, one grid is enough
        bitset_grid rw_grid_;
        // Location we are evaluating agents on
        location focus_;
        // size of bitsets
        uint size_;
        // even / uneven partitioning of grid for diffusion
        uint period_;
        // the original bitsets
        std::vector< boost::dynamic_bitset<> > bits_;
        
        protected:
        static double diffusion_;
    };
    
    inline void BitEnvironment::getBits( 
        boost::dynamic_bitset<> &b, const location &loc ) const
    { b = rw_grid_[ loc.first ][ loc.second ]; }

    inline void BitEnvironment::setBits(
        const boost::dynamic_bitset<> &b, const location &loc ) 
    { rw_grid_[ loc.first ][ loc.second ] = b; }

    inline const BitEnvironment::bitset_grid &
    BitEnvironment::grid() const
    { return rw_grid_; }
    
    inline const std::vector< boost::dynamic_bitset<> > &
    BitEnvironment::bits() const
    { return bits_; }
    
    /// \class NoisyBitEnvironment
    /// \brief Environment with bitstrings and noise on them
    class NoisyBitEnvironment : public BitEnvironment {
        public:
        // Constructor
        NoisyBitEnvironment();
        // Constructor w/ dimensions
        NoisyBitEnvironment( int, int, uint, double );
        // Destructor
        ~NoisyBitEnvironment() {};

        /// Initialise grid with a few bitsets
        virtual void initialise( 
            const std::vector< boost::dynamic_bitset<> > & );
        /// Fluctuations are noise
        virtual void fluctuate( long );
        /// The environment is changed by an agent
        virtual void update( const Agent *, const location & );
        
        private:
        void generateLocations();
        
        private:
        double rate_;
        std::vector< location > cache_locs_;
    };

    
    /// \class NoRotateBitEnvironment
    /// \brief Environment with bitstrings that can be exhausted
    class NoRotateBitEnvironment : public BitEnvironment {
        public:
        // Constructor
        NoRotateBitEnvironment();
        // Constructor w/ dimensions
        NoRotateBitEnvironment( int, int, uint, double );
        // Destructor
        ~NoRotateBitEnvironment() {};

        /// Initialise grid with bitset and intro heterogeneity
        virtual void initialise( 
            const std::vector< boost::dynamic_bitset<> > & );
        /// Fluctuations are influx of new bitstrings (only at empty sites)
        virtual void fluctuate( long );
        /// The environment is changed by an agent
        virtual void update( const Agent *, const location & );
        
        private:
        double influx_;
    };
    
    /// \class SubsetBitEnvironment
    /// \brief Environment with bitstrings that are not changed, only read
    class SubsetBitEnvironment : public BitEnvironment {
        public:
        // Constructor
        SubsetBitEnvironment();
        // Constructor w/ dimensions
        SubsetBitEnvironment( int, int, uint );
        // Destructor
        ~SubsetBitEnvironment() {};

        /// Initialise grid with bitset and intro heterogeneity
        virtual void initialise( 
            const std::vector< boost::dynamic_bitset<> > & );
        /// Fluctuations are non-existent
        virtual void fluctuate( long );
        /// The environment is changed by an agent
        virtual void update( const Agent *, const location & ) {};
    };
    
    /// \class OneBiteBitEnvironment
    /// \brief Environment with bitstrings that are not changed, only read
    class OneBiteBitEnvironment : public BitEnvironment {
        public:
        // Constructor
        OneBiteBitEnvironment();
        // Constructor w/ dimensions
        OneBiteBitEnvironment( int, int, uint, double );
        // Destructor
        ~OneBiteBitEnvironment() {};

        /// Initialise grid with bitset and intro heterogeneity
        virtual void initialise( 
            const std::vector< boost::dynamic_bitset<> > & );
        /// Fluctuations are non-existent
        virtual void fluctuate( long );
        /// The environment is changed by an agent
        virtual void update( const Agent *, const location & );
        
        private:
        double influx_;
    };
}
#endif
