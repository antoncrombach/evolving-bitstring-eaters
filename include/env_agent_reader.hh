//
// XML SAX reader for agents
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_ENV_AGENT_READER_
#define _DORSALFIN_ENV_AGENT_READER_

#include "defs.hh"
#include "env_agent.hh"
#include "gene.hh"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace dorsalfin {

// Xerces shit.
XERCES_CPP_NAMESPACE_USE

XERCES_CPP_NAMESPACE_BEGIN
    class AttributeList;
XERCES_CPP_NAMESPACE_END

    /// \class EnvAgentReader
    /// \brief Read in an agent from an xml and dot file.
    ///
    /// EnvAgentReader reads a genome from an xml file and a network from a dot
    /// file.
    class EnvAgentReader : public DefaultHandler {
        public:
        /// Constructor.
        EnvAgentReader( Factory *, int );
        /// Destructor.
        ~EnvAgentReader() {};

        // Return the agent
        GenRegAgent* agent();
        
        /// Read opening xml tag.
        void startElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname,
            const Attributes &attrs );
        /// Read tag contents. 
        void characters( const XMLCh* const chars, const uint length );
        /// Read closing xml tag.
        void endElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname );
        /// Set stuff to zero 
        void resetDocument();
    
        /// A few functions for catching warnings and errors
        void error( const SAXParseException & );
        void fatalError( const SAXParseException & );
        void resetErrors();
        bool sawErrors() const;
        bool done() const;
        
        private:
        Factory *factory_;
        bool sawErrors_, done_, class_;
        
        std::string agent_class_;
        int type_;
        
        GenRegAgent *agent_;
        AgentTag me_, mother_;
        Network *network_;
        Genome *genome_;
        Chromosome *chromo_;
        std::list< ChromosomeElement* > *chr_;
        std::set< int > gene_set_;
        boost::dynamic_bitset<> reference_;
    };

inline bool EnvAgentReader::sawErrors() const
{ return sawErrors_; }

inline bool EnvAgentReader::done() const
{ return done_; }

inline void EnvAgentReader::error( const SAXParseException &e ) 
{ sawErrors_ = true; }

inline void EnvAgentReader::fatalError( const SAXParseException &e )
{ sawErrors_ = true; }

inline void EnvAgentReader::resetErrors()
{ sawErrors_ = false; }

inline GenRegAgent* EnvAgentReader::agent() 
{ return agent_; }
}
#endif
