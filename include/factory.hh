//
// Factory for creating a genome and other model entities.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_FACTORY_H_
#define _DORSALFIN_FACTORY_H_

#include "defs.hh"
#include "config.hh"

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/XMLString.hpp>


namespace dorsalfin {

    XERCES_CPP_NAMESPACE_USE
    
    /// \class Factory
    /// \brief Model parts are built in the \c Factory.
    ///
    /// Given the configuration, parts of the model can be built in the factory.
    class Factory {
        
        public:
        /// (Dummy) constructor
        Factory();
        /// Constructor with model and configuration
        Factory( Config *g );
        /// Destructor
        ~Factory();

        /// Create a gene
        Gene* gene();
        /*/// Create a transcription factor
        TransFactor* transfactor();*/
        /// Create a transposon (a special type of genes).
        Retroposon* retroposon();
        /// Create a repeat (long terminal repeat)
        Repeat* repeat();
        
        /// Create a chromosome. It is filled with binding sites,
        /// genes and transposons (not interrupting any upstream)
        /// according to the settings in the configuration
        Chromosome* chromosome( int );
        /// Create a genome
        Genome* genome( int );
        /// Provide network w/ information (such as reference network)
        void network( int );

        /// Create an agent. The type of agent is given in the 
        /// configuration
        Agent* agent( int );
        /// Create the population, still needs to be initialised.
        Population* population();
        /// Create the environment, still needs to be initialised.
        /// The type is given in the configuration.
        Environment* environment();

        /// Create a scaling scheme used to scale raw scores of the agents
        ScalingScheme* scalingScheme();
        /// Create a selection scheme used to select an agent for
        /// reproduction based on its score.
        SelectionScheme* selectionScheme();

        /// Create an observer manager. It keeps track of observers and to
        /// which subject they are linked.
        ObserverManager* observerManager();

        /// Create a queue of events
        std::priority_queue< Event* >* eventQueue();

        /// Configuration used for creating agents etc
        const Config & configuration() const;
        
        private:
        // Read an agent from file (xml)
        GenRegAgent* readEnvAgent( const std::string &, int );
        // Read an entire population from file (xml)
        boost::tuple< std::vector< Agent* >, std::vector< location > > 
            readEnvPop( const std::string &, int );
        // Read environment
        Environment* readEnvironment( const std::string &,
            std::vector< boost::dynamic_bitset<> > & );
        
        // Read the bits of environmental attractors
        void readInput( std::vector< boost::dynamic_bitset<> > & );
        // Read the agent's config file 
        void readAgentConfigurations();
        // Set some configuration things
        void setConfiguration( int );
        
        // Build correct # of elements according to config
        std::vector< std::list< ChromosomeElement* > > *
            chromosomeParts( int );
        // Organise chromosome up to certain degree
        std::list< ChromosomeElement* >* organiseChromosome( 
            std::vector< std::list< ChromosomeElement* > > *, double );
        
        // Reset the tag counter
        void resetTag();
        // Produce a new id
        int tag();

        private:
        Config *conf_;
        int gene_tag_;
    };

    inline const Config & Factory::configuration() const
    { return *conf_; }
    
    inline void Factory::resetTag()
    { gene_tag_ = 0; }
    
    inline int Factory::tag()
    { return gene_tag_++; }
}
#endif

