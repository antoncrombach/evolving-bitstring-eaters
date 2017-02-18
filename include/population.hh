//
// Interface of a population (may be refactored to an abstract class etc).
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_POPULATION_H_
#define _DORSALFIN_POPULATION_H_

#include "defs.hh"
#include "observer.hh"
#include "subject.hh"
#include "model.hh"
#include "agent.hh"
#include "scaling.hh"
#include "selection.hh"

namespace dorsalfin {

    /// \class Population
    /// \brief Collection of agents on a grid.
    ///
    /// Agents live on a grid. They reproduce into empty grid locations, 
    /// which depends on the fitness score. The exact nature of the fitness
    /// is given by the agent.
    class Population : public Subject, public XmlWriter {
        public:
        /// Typedef of the agent grid
        typedef boost::multi_array< Agent*, 2 > agents_grid;
        /// Mapping agents on the grid (inverse of a grid)
        typedef std::map< Agent*, location > agents_map;
        /// Range of agent grid
        typedef agents_grid::index_range range;
        /// Agents map iterator
        typedef agents_map::iterator map_ag_iter;
        /// Const agents map iterator
        typedef agents_map::const_iterator const_map_ag_iter;
        /// Vector of agents iterator
        typedef std::vector< Agent* >::iterator ag_iter;
        /// Vector of locations iterator
        typedef std::vector< location >::iterator loc_iter;
        /// Const vector of locations iterator
        typedef std::vector< location >::const_iterator const_loc_iter;
        
        public:
        /// Enumeration for synchronous updating
        enum grid_type { reading, writing };
        
        public:
        /// (Dummy) constructor
        Population();
        /// Constructor with grid dimensions, initial agents, a scaling
        /// scheme and a selection scheme.
        Population( int, int, std::vector< Agent* > &, 
            ScalingScheme*, SelectionScheme* );
        /// Constructor with locations for the agents
        Population( int, int, std::vector< Agent* > &, 
            const std::vector< location > &, ScalingScheme*,
            SelectionScheme* );
        /// Copy constructor (deep copy)
        Population( const Population & );
        /// Destructor
        ~Population();
        
        /// Get the moore neighbourhood of the given location (from shadow)
        std::vector< Agent* > moore( const location & );
        /// Get the moore neighbourhood of the agent (from shadow)
        std::vector< Agent* > moore( const Agent & );
        /// Get the location of the agent (from shadow)
        location whereIs( const Agent & ); 
        
        /// Tag all the agents such that we can identify them
        void initialise();
        /// Evaluate all agents (in shadow) on the environment
        void evaluate( const Environment & );
        /// Perform one update, one timestep
        void step();
        /// Do some things to end the simulation nicely
        //void finish();

        /// Insert several agents at random locations (in grid)
        void insert( std::vector< Agent* > & );
        /// Insert several agents at specified locations (in grid)
        void insertAt( ag_iter, ag_iter, const std::vector< location > & );
        /// Insert several agents at specified locations (in grid)
        void insertAt( std::vector< Agent* > &, 
            const std::vector< location > & );
        /// Insert an agent at a location (in grid)
        void insertAt( Agent*, const location & );
        /// Erase the given agents (from grid)
        void erase( std::vector< Agent* > & );
        /// Erase the agent (from grid)
        void erase( Agent* );
        /// Erase any agent at the given locations (from grid)
        void eraseAt( const std::vector< location > & );
        /// Erase any agent at the specified location (from grid)
        void eraseAt( location );
        
        /// Get the grid. Only used by observers
        const agents_grid & grid() const;
        /// Get the map
        const agents_map & map() const;
        /// Get the current generation number (time)
        long generation() const;
        /// Set the model this population belongs to
        void model( Model * );
        /// Get the model this population belongs to
        Model* model() const;
        
        /// Get the current number of agents 
        uint nrAgents() const;
        /// Is the grid empty?
        bool empty() const;
        /// Special function to see if only one type of agents is left
        bool hasEveryAgentType() const;
        
        /// Overload a subject method
        void attach1( AsyncLogObserver * );
        void attach2( AsyncLogObserver * );
        void attach3( AsyncLogObserver * );
        /// And another method for observers
        void closeAll();
        /// Breaching the idea of not knowing the exact nature of observers...
        void initAncestorTracing();
        
        /// write a nice xml formatted version
        virtual void writeXml( std::ostream & ) const;
        
        public:
        /// Set the threshold for not performing an action due to too few
        /// agents in the neighbourhood. (Already scaled value should be given.)
        static void threshold( double );
        /// Get the threshold for not performing an action due to too few
        /// agents in the neighbourhood
        static double threshold();
        /// Set the number of different agents on the plane
        static void nrAgentTypes( int );
        /// Get the number of agent types
        static int nrAgentTypes();
        /// Set the placement type (random, patch)
        static void placement( const std::string & );
        /// Get the placement
        static std::string placement();
        /// Set if we shuffle or not
        static void shuffling( bool );
        /// Are we shuffling?
        static bool shuffling();
        
        private:
        class IsNotAvailable :
            public std::unary_function< Agent*, bool > {
            public:
                IsNotAvailable() {}
                bool operator()( Agent *ce ) const 
                { return ce == static_cast< Agent* >( 0 ); }
        };
        
        struct SmallerTypeThan :
            public std::binary_function< Agent*, Agent*, bool > {
                bool operator()( Agent *ag1, Agent *ag2 ) const 
                { return ag1->type() < ag2->type(); }
        };
        
        private:
        // get some random empty locations
        std::vector< location > randEmptyLocations( int );
        // get a patch (squarish) of empty locations
        std::vector< std::vector< location > > 
            patchEmptyLocations( std::vector< int > );
        // auxilary function filling a vector with all locations in the grid
        void asLocations( std::vector< location > & );
        // which plane to zero?
        void zero( grid_type );
        // swap the two planes and make them consistent
        void swap();
        // swap two locations on a plane
        void swap( grid_type, const location &, const location & );
        // shuffle the locations of both planes
        void shuffle();
        // do not erase agent, only its pointers in the write plane and map
        void shallowErase( Agent * );
        // do not erase agent, only its pointers...
        void shallowEraseAt( location );
        
        private:
        // read from shadow_grid, write to grid
        agents_grid plane_one_, plane_two_;
        agents_grid *write_grid_, *read_grid_;
        agents_map write_agents_, read_agents_;
        
        std::vector< location > shuffle_locs_;
        
        ScalingScheme *scaling_;
        SelectionScheme *selection_;
        
        Model *model_;
        AsyncLogObserver *async_agent_obs_;
        AsyncLogObserver *async_env_change_;
        // Says DSBs but logs all kinds of mutations
        AsyncLogObserver *async_dsbs_;
        
        private:
        static bool shuffle_;
        static double threshold_;
        static int nr_agent_types_;
        static std::string placement_;
        
        const int agent_view_size_;
    };

    
    inline const Population::agents_grid & Population::grid() const
    { return *write_grid_; }

    inline const Population::agents_map & Population::map() const
    { return write_agents_; }

    inline long Population::generation() const
    { return model_->now(); }

    inline void Population::model( Model *m )
    { model_ = m; }

    inline Model* Population::model() const
    { return model_; }
    
    inline uint Population::nrAgents() const
    { return write_agents_.size(); }
    
    inline bool Population::empty() const
    { return write_agents_.empty(); }
}
#endif


