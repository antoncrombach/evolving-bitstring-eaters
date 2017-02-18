//
// XML SAX reader for populations
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_ENV_POP_READER_
#define _DORSALFIN_ENV_POP_READER_

#include "defs.hh"
#include "env_agent.hh"
#include "gene.hh"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace dorsalfin {

XERCES_CPP_NAMESPACE_USE

XERCES_CPP_NAMESPACE_BEGIN
    class AttributeList;
XERCES_CPP_NAMESPACE_END

    /// \class EnvPopReader
    /// \brief Read a population of agents from an xml file and a set of dot
    /// files.
    class EnvPopReader : public DefaultHandler {
        public:
        /// ctor
        EnvPopReader( Factory *, int );
        /// dtor
        ~EnvPopReader() {};

        /// Return the population
        boost::tuple< std::vector< Agent* >, std::vector< location > >
            population();
        
        /// Read opening tag
        void startElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname,
            const Attributes &attrs );
        /// Read tag contents (not used a lot)
        void characters( const XMLCh* const chars, 
            const unsigned int length );
        /// Read closing tag
        void endElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname );
        /// Stuff to zero
        void resetDocument();
    
        /// A few functions for catching warnings and errors
        void error( const SAXParseException &exc );
        void fatalError( const SAXParseException &exc );
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
        Genome *genome_;
        Chromosome *chromo_;
        std::list< ChromosomeElement* > *chr_;
        
        std::vector< Agent* > population_;
        std::vector< location > locations_;
        std::set< int > gene_set_;
        boost::dynamic_bitset<> reference_;
    };

inline bool EnvPopReader::sawErrors() const
{ return sawErrors_; }

inline bool EnvPopReader::done() const
{ return done_; }

inline void EnvPopReader::error( const SAXParseException &e ) 
{ sawErrors_ = true; }

inline void EnvPopReader::fatalError( const SAXParseException &e )
{ sawErrors_ = true; }

inline void EnvPopReader::resetErrors()
{ sawErrors_ = false; }

inline
boost::tuple< std::vector< dorsalfin::Agent* >, std::vector< dorsalfin::location > >
EnvPopReader::population() 
{ return boost::make_tuple( population_, locations_ ); }
}
#endif
