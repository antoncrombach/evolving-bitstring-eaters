//
// Hopfield directed graph
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_HOPFIELD_GRAPH_H_
#define _DORSALFIN_HOPFIELD_GRAPH_H_

#include "defs.hh"
#include "ref_node.hh"

namespace dorsalfin {

    /// \class Network
    /// \brief Gene regulation network
    ///
    /// Important: at reading in the starting (reference) network, any gene
    /// labelled 'transfac' is a transcription factor during the simulation.
    class Network : public XmlWriter {
        public:
        /// Typedef edges
        typedef std::pair< int, int > dorsalfin_edge;
        /// Typedef of the graph storage structure
        typedef boost::adjacency_list< 
            boost::vecS, boost::vecS, boost::bidirectionalS, 
            boost::property< boost::vertex_name_t, Gene* >,
            boost::property< boost::edge_weight_t, int > > dorsalfin_graph;
        typedef boost::adjacency_list< 
            boost::vecS, boost::vecS, boost::bidirectionalS, RefNode, 
            boost::property< boost::edge_weight_t, int > > ref_graph;
        /// Property maps
        /// Network/graph
        typedef boost::property_map< dorsalfin_graph, 
            boost::vertex_name_t >::type name_map;
        typedef boost::property_map< dorsalfin_graph, 
            boost::edge_weight_t >::type interaction_map;
        /// Accessing vertices and edges
        typedef boost::graph_traits< dorsalfin_graph >::vertex_iterator
            node_iter;
        typedef boost::graph_traits< dorsalfin_graph >::in_edge_iterator
            edge_iter;
        /// and for the reference graph
        typedef boost::graph_traits< ref_graph >::vertex_iterator ref_node_iter;
        typedef boost::graph_traits< ref_graph >::edge_iterator ref_edge_iter;
        // and for the edge mapping
        typedef std::map< dorsalfin_edge, int > edge_iaction_map;
        typedef edge_iaction_map::iterator edge_iaction_iter;
        typedef edge_iaction_map::const_iterator edge_iaction_citer;
        
        public:
        /// Standard ctor
        /// pre: reference network should be read in already
        Network();
        /// Copy ctor
        explicit Network( const Network & );
        /// Cloning
        Network* clone() const;
        /// Copying method
        void copy( const Network & );
        /// Dtor
        ~Network();

        /// Build the graph from the genome
        void build( const Genome & );
        /// Propagate the network, stop if output mismatches input
        void propagate();
        /// Hamming distance to other state
        int hamming( const boost::dynamic_bitset<> & ) const;
        /// Hamming distance of the path of the output node to the other bitset
        int pathHamming( const boost::dynamic_bitset<> & ) const;
        
        /// Do we have an output node?
        bool hasOutputNode() const;
        /*/// Are we in an attractor?
        bool inAttractor() const;*/
        /// Are all nodes zero?
        bool allZero() const;
        /// Is network already build?
        bool isEmpty() const;
        
        /// Get graph
        const dorsalfin_graph & graph() const;
        /// Get output
        boost::dynamic_bitset<> output() const;
        /// State as a bitset
        void state( boost::dynamic_bitset<> & ) const;
        
        /// Get the input bitset
        boost::dynamic_bitset<> input() const;
        /// Set all genes to zero
        void resetState();
        /// Set the input bitset
        void setInput( const boost::dynamic_bitset<> & );
        /// Set the state of input genes
        void setInputState( const boost::dynamic_bitset<> & );
        
        /// Write in dot format (although it is called xml)
        virtual void writeXml( std::ostream & ) const;

        public:
        /// Set lead-in time
        static void initPropagate( int );        
        /// Set time limit
        static void maxPropagate( int );
        /// Get time limit
        static int maxPropagate();
        /// Sequential or parallel update
        static void seqPropagate( bool );
        /// Write reference graph to file
        static void writeReferenceGraph( std::ostream & );
        /// Read reference graph from file
        static void readReferenceGraph( std::istream & );
        /// Get nr of input genes (assuming consecutive numbers from zero)
        static int referenceNrInputIds();
        /// Get reference node ids of genes
        static const std::set< int > & referenceIds();
        /// Get reference node ids of input genes
        static const std::set< int > & referenceInputIds();
        /// Get reference node ids of transcription factors
        static const std::set< int > & referenceTransFacIds();
        /// Get reference nodes
        static const std::map< int, RefNode > & referenceNodes();
        /// Get reference edges
        static const std::map< dorsalfin_edge, int > & referenceEdges();
        
        private:
        // Set dynamic_bitset \c state_ to current network state
        void doState();
        // Set state of output (pushed back to output_path_)
        bool outputState();

        private:
        // Auxilary function for reading the reference network from file
        static void doReadReferenceGraph( std::istream &is, ref_graph &ref,
            boost::dynamic_properties &dp );
        
        protected:
        // network of genes
        mutable dorsalfin_graph graph_;
        // state of network and output produced by it
        boost::dynamic_bitset<> state_, output_path_, input_path_;
        // caching 
        name_map graph_gmap_;
        interaction_map graph_imap_;
        // signal if output node present, if in attractor
        bool has_output_, in_attractor_;
        
        protected:
        // Lead-in time for network calculation
        static int init_propagate_;
        // Time limit for stable state calculation
        static int max_propagate_;
        // Sequential update (true) or parallel (false)
        static bool seq_propagate_;
        // How many nodes do we have?
        static int nr_nodes_;
        // Essential nodes (used for checking if lethal, ie missing a node)
        // Right now, it is only the output node
        static std::set< int > ref_output_ids_;
        // Nodes that are input nodes
        static std::set< int > ref_input_ids_;
        // Nodes that are transcription factors (input nodes incl)
        static std::set< int > ref_transfac_ids_;
        // All nodes
        static std::set< int > ref_ids_;
        // All nodes with their info
        static std::map< int, RefNode > ref_nodes_;
        // Which edges have we got (structure of net is predefined)
        static std::map< dorsalfin_edge, int > ref_edges_;
    };
    
    inline bool Network::hasOutputNode() const
    { return has_output_; }
    
    /*inline bool Network::inAttractor() const
    { return in_attractor_; }*/
    
    inline bool Network::isEmpty() const
    { return boost::num_vertices( graph_ ) == 0; }
    
    inline const Network::dorsalfin_graph & Network::graph() const
    { return graph_; }
    
    inline void Network::setInput( const boost::dynamic_bitset<> &b )
    { input_path_ = b; }
    
    inline boost::dynamic_bitset<> Network::input() const
    { return input_path_; }
    
    inline boost::dynamic_bitset<> Network::output() const
    { return output_path_; }
    
    inline void Network::state( boost::dynamic_bitset<> &v ) const 
    { v = state_; }
}
#endif
