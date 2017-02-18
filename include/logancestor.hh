//
// Some observers of the population.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_LOGANCESTOR_H_
#define _DORSALFIN_LOGANCESTOR_H_

#include "defs.hh"
#include "observer.hh"
#include "population.hh"
#include "stream_manager.hh"

namespace dorsalfin {

    /// \class LogCsvAncestors
    ///
    /// Log ancestors to trace them back
    class LogCsvAncestors : public AsyncLogObserver {
        public:
        // ancestor graph
        typedef boost::adjacency_list< 
            boost::setS, boost::setS, boost::bidirectionalS,
            boost::property< boost::vertex_index_t, uint, AgentTag >,
            boost::property< boost::edge_index_t, uint > > ancestor_graph;
        // vertex descriptor
        typedef ancestor_graph::vertex_descriptor anc_vertex;
        typedef ancestor_graph::edge_descriptor anc_edge;
        // a mapping
        typedef boost::property_map< ancestor_graph, 
            boost::vertex_bundle_t >::type anc_birth_map;
        // accessing vertices
        typedef boost::graph_traits< ancestor_graph >::vertex_iterator
            ver_iter;
        // accessing edges
        typedef boost::graph_traits< ancestor_graph >::adjacency_iterator
            adj_iter;
        typedef boost::inv_adjacency_iterator_generator< ancestor_graph >::type
            invadj_iter;
            
        public:
        LogCsvAncestors( std::string );
        virtual ~LogCsvAncestors() {}
        
        virtual void initialize( Subject * );
        virtual void update( Subject * );
        virtual void finalize();
        
        public:
        static void maxForestSize( uint );
        
        private:
        void writeHeader();
        void writeForest();
        
        void renewForest();
        void initForest();
        
        void birth( Agent * );
        void death( Agent * );
        
        void mapToTime( std::map< AgentTag, anc_vertex > & );
        
        private:
        int nr_roots_;
        anc_vertex root_;
        ancestor_graph forest_;
        std::map< AgentTag, anc_vertex > leafs_, mothers_;
        uint vertex_nr_, edge_nr_;
        
        private:
        static uint max_forest_size_;
    };

    /// \class LogXmlAgentTrace
    ///
    /// Log agents given a trace of agents to look for...
    class LogXmlAgentTrace : public AsyncLogObserver {
        // The trace file has a certain format. Unfortunately xml is a bit of
        // a deception, so the trace file is in csv format:
        //
        // each line is an agent id-tag: birth <tab> x <tab> y <nl>
        public:
        LogXmlAgentTrace( std::string, std::string );
        virtual ~LogXmlAgentTrace() {}
        
        virtual void initialize( Subject * );
        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        void writeFooter();
        void initNextTarget();
        void nextTargets();
        
        private:
        bool first_;
        std::string dname_;
        AgentTag next_target_;
        std::set< AgentTag > targets_;
        boost::filesystem::ifstream *trace_;
    };
}
#endif

