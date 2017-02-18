//
// XML sax reader for populations
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "defs.hh"
#include "factory.hh"

#include "env_pop_reader.hh"

#include "population.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "retroposon.hh"
#include "repeat.hh"

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

dorsalfin::EnvPopReader::EnvPopReader( Factory *ft, int type ) 
    : factory_( ft ), 
      sawErrors_( false ), done_( false ), class_( false ),
      agent_class_(), type_( 0 ), 
      agent_( 0 ), genome_( 0 ), chromo_( 0 ), chr_( 0 ),
      population_( std::vector< Agent* >() ), 
      locations_( std::vector< location >() ), gene_set_() {}

void
dorsalfin::EnvPopReader::startElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname,
        const Attributes &attrs ) {
    // now check type of element
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // create tag
        char *bux = XMLString::transcode( attrs.getValue( 0u ) );
        me_.time = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        me_.loc.first = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        me_.loc.second = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 3 ) );
        me_.i = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // and create a new chromosome list
        chr_ = new std::list< ChromosomeElement* >();
    } else if( aux == "class" ) {
        // we know what kind of agent to create
        class_ = true;
    } else if( aux == "parent" ) {
        // create parent tag
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        mother_.time = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        mother_.loc.first = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        mother_.loc.second = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );        
        bux = XMLString::transcode( attrs.getValue( 3 ) );
        mother_.i = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
    } else if( aux == "gene" ) {
        // create another gene
        // id
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // state
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        int dux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // threshold
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        int eux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // create gene
        chr_->push_back( new Gene( cux, dux, eux ) );
        // and count it
        gene_set_.insert( cux );
    } else if( aux == "repeat" ) {
        // create repeat
        chr_->push_back( new Repeat() );
    } else if( aux == "tposon" ) {
        // create tposon
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        chr_->push_back( new Retroposon( cux ) );
    } else if( aux == "ia" ) {
        // create interaction
        // target
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // interaction value (weight)
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        int dux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );        
        chr_->push_back( new Interaction( cux, dux ) );
    }
}

void
dorsalfin::EnvPopReader::characters( const XMLCh* const chars, 
        const unsigned int length) {
    char *xxx = XMLString::transcode( chars );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( class_ ) {
        agent_class_ = aux;
        class_ = false;
    }
}

void
dorsalfin::EnvPopReader::endElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname ) {
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // ready to create the right agent
        if( agent_class_ == "EnvAgent" ) {
            agent_ = new EnvAgent( type_, genome_ );
            agent_->myTag( me_ );
            agent_->parentTag( mother_ );
            agent_->readFromFile( true );
            // and insert in population
            population_.push_back( agent_ );
            locations_.push_back( me_.loc );
        } else {
            throw "Don't know this agent! Peace out..";
        }
        resetDocument();
    } else if( aux == "genome" ) {
        // create genome (with only one chromosome)
        std::list< Chromosome* > *ll = new std::list< Chromosome* >();
        ll->push_back( chromo_ );
        genome_ = new Genome( ll );
    } else if( aux == "chromosome" ) {
        // create chromosome
        const Config &cf = factory_->configuration();
        chromo_ = new Chromosome( 0, chr_ );
        chromo_->copyGeneRate( cf.optionAsDouble( "cp_gene", type_ ) );
        chromo_->removeGeneRate( cf.optionAsDouble( "rm_gene", type_ ) );
        chromo_->thresholdGeneRate( cf.optionAsDouble( "thr_gene", type_ ) );
        chromo_->tagGeneRate( cf.optionAsDouble( "tag_gene", type_ ) );
        chromo_->stateRate( cf.optionAsDouble( "state_gene", type_ ) );
        
        chromo_->copyRetroposonRate( cf.optionAsDouble( "cp_tp", type_ ) );
        chromo_->removeRetroposonRate( cf.optionAsDouble( "rm_tp", type_ ) );
        chromo_->newRetroposonRate( cf.optionAsDouble( "new_tp", type_ ) );
        chromo_->removeRepeatRate( cf.optionAsDouble( "rm_ltr", type_ ) );
        chromo_->recombinationRate( 
            cf.optionAsDouble( "dsb_recombination", type_ ) );
        
        chromo_->weightInteractionRate( 
            cf.optionAsDouble( "weight_ia", type_ ) );
        chromo_->copyInteractionRate( cf.optionAsDouble( "cp_ia", type_ ) );
        chromo_->removeInteractionRate( cf.optionAsDouble( "rm_ia", type_ ) );
        chromo_->newInteractionRate( cf.optionAsDouble( "new_ia", type_ ) );
        chromo_->tagInteractionRate( cf.optionAsDouble( "tag_ia", type_ ) );
    } else if( aux == "simulation" ) {
        done_ = true;
    }
}

void
dorsalfin::EnvPopReader::resetDocument() {
    // just set to zero, we do not own this population
    agent_ = 0;
    genome_ = 0;
    chromo_ = 0;
    chr_ = 0;
    class_ = false;
    gene_set_.clear();
    reference_.clear();
}

