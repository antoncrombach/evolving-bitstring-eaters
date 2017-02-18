//
// Implementation of the model factory
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "factory.hh"
#include "stream_manager.hh"
#include "event.hh"

#include "environment.hh"
#include "population.hh"
#include "simple_agent.hh"
#include "env_agent.hh"

#include "genome.hh"
#include "chromosome.hh"
#include "gene.hh"
#include "tfactor.hh"
#include "interaction.hh"
#include "retroposon.hh"
#include "repeat.hh"

#include "network.hh"

#include "scaling.hh"
#include "selection.hh"

#include "env_agent_reader.hh"
#include "env_pop_reader.hh"

XERCES_CPP_NAMESPACE_USE

dorsalfin::Factory::Factory() {
    gene_tag_ = 0;
}

dorsalfin::Factory::Factory( Config *gc ) : conf_( gc ) {
    gene_tag_ = 0;
}

dorsalfin::Factory::~Factory() {}

dorsalfin::Gene*
dorsalfin::Factory::gene() {
    return new Gene( tag() );
}

/*dorsalfin::TransFactor*
dorsalfin::Factory::transfactor() {
    return new TransFactor( tag() );
}*/

dorsalfin::Retroposon* 
dorsalfin::Factory::retroposon() {
    return new Retroposon( tag() );
}

dorsalfin::Repeat*
dorsalfin::Factory::repeat() {
    return new Repeat();
}

dorsalfin::Chromosome* 
dorsalfin::Factory::chromosome( int type ) {
    // auxilary typedef
    std::list< ChromosomeElement* > *ll = 
        organiseChromosome( chromosomeParts( type ), 
            conf_->optionAsDouble( "organised", type ) );
    
    // scan list and flank all retroposons with repeats
    Chromosome::ce_iter i = ll->begin();
    while( i != ll->end() ) {
        if( Chromosome::IsRetroposon()( *i ) ) {
            ll->insert( i, repeat() );
            ll->insert( boost::next( i ), repeat() );
        }
        ++i;
    }
    
    // init with right values for this type of agent
    Chromosome *result = new Chromosome( 0, ll );
    return result;
}

dorsalfin::Genome* 
dorsalfin::Factory::genome( int type ) {
    resetTag();
    // #chromosomes
    int n_ch = conf_->optionAsInt( "chromos", type );

    std::list< Chromosome* > *ll = new std::list< Chromosome* >();
    for( int i = 0; i < n_ch; ++i ) {
        ll->push_back( this->chromosome( type ) );
    }
    return new Genome( ll );   
}

void
dorsalfin::Factory::network( int type ) {
#ifdef DEBUG
    cout << "! init reference network" << endl;
#endif
    // read reference network
    boost::filesystem::ifstream *file = 
        StreamManager::instance()->openInFileStream( conf_->optionAsString( "network", type ), boost::filesystem::fstream::in );
    Network::readReferenceGraph( *file );
    StreamManager::instance()->closeInFileStream( file );
    Network::initPropagate( conf_->optionAsInt( "init_propagate", type ) );
    Network::maxPropagate( conf_->optionAsInt( "max_propagate", type ) );

    // write reference network
    boost::filesystem::ofstream *outfile = StreamManager::instance()->
        openOutFileStream( std::string( "refnet.dot" ),
        boost::filesystem::fstream::out );
    Network::writeReferenceGraph( *outfile );
    StreamManager::instance()->closeOutFileStream( outfile );
}

dorsalfin::Agent*
dorsalfin::Factory::agent( int type ) {
    // build agent according to type    
    Agent *result;
    if( conf_->optionAsString( "agent", type ) == "simple" ) {
        result = new SimpleAgent();
    } else if( conf_->optionAsString( "agent", type ) == "env" ) {
        result = new EnvAgent( type, genome( type ) );
        network( type );
    } else {
        // try to open it as a file
        result = readEnvAgent( conf_->optionAsString( "agent", type ), type );
        network( type );
    }
    // assert
    if( static_cast< int >( Network::referenceIds().size() ) != 
        conf_->optionAsInt( "genes", type ) ) {
        std::cout << "Warning: network size != nr of genes" << std::endl;
    }
    // and continue
    result->initialise();
    return result;
}

dorsalfin::Population*
dorsalfin::Factory::population() {
    Population *result = 0;
    // read some population level parameters
    Population::nrAgentTypes( conf_->optionAsInt( "nr_agent_type" ) );
    Population::placement( conf_->optionAsString( "agent_placement" ) );
    Population::shuffling( conf_->optionAsString( "shuffle" ) == "true" );
    // and per agent type stuff
    readAgentConfigurations();
    
    uint xx = conf_->optionAsInt( "grid_x" );
    uint yy = conf_->optionAsInt( "grid_y" );
    bool first_pop = conf_->hasOption( "population_one" );
    bool second_pop = conf_->hasOption( "population_two" );
    bool rand_pop = conf_->optionAsString( "population_start" ) == "random"; 
    if( !first_pop and !second_pop ) {
        std::cout << "Creating a population" << std::endl;
        std::vector< Agent* > aux;
        int cux = conf_->optionAsInt( "nr_agent_type" );
        // random start or homogeneous start
        if( rand_pop ) {
            for( int i = 0; i != cux; ++i ) {
                int bux = conf_->optionAsInt( "init_nr_agents", i );
                for( int j = 0; j != bux; ++j ) {
                    aux.push_back( agent( i ) );
                }
            }
        } else {
            for( int i = 0; i != cux; ++i ) {
                int bux = conf_->optionAsInt( "init_nr_agents", i );
                if( bux > 0 ) {
                    aux.push_back( agent( i ) );
                    for( int j = 1; j != bux; ++j ) {
                        aux.push_back( aux.front()->clone() );
                    }
                }
            }
        }
        result = new Population( xx, yy, aux,
            scalingScheme(), selectionScheme() );
    } else if( first_pop and !second_pop ) {
        std::cout << "Reading in one population (1)" << std::endl;
        std::vector< Agent* > aux;
        std::vector< location > loc;
        boost::tie( aux, loc ) = 
            readEnvPop( conf_->optionAsString( "population_one" ), 0 );
        // set agent type to 0 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            //( **ii ).type( 0 );
            ( **ii ).initialise();
        } 
        result = new Population( xx, yy, aux, loc,
            scalingScheme(), selectionScheme() );
    } else if( !first_pop and second_pop ) {
        std::cout << "Reading in one population (2)" << std::endl;
        std::vector< Agent* > aux;
        std::vector< location > loc;
        boost::tie( aux, loc ) = 
            readEnvPop( conf_->optionAsString( "population_two" ), 0 );
        // set agent type to 0 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            //( **ii ).type( 0 );
            ( **ii ).initialise();
        } 
        result = new Population( xx, yy, aux, loc,
            scalingScheme(), selectionScheme() );
    } else { // both populations
        std::cout << "Reading in two populations" << std::endl;
        std::vector< Agent* > aux;
        std::vector< location > loc;
        boost::tie( aux, loc ) = 
            readEnvPop( conf_->optionAsString( "population_one" ), 0 );
        // set agent type to 0 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            ( **ii ).type( 0 );
            ( **ii ).initialise();
        } 
        // 2nd population
        std::vector< Agent* > bux;
        std::vector< location > bloc;
        boost::tie( bux, bloc ) = 
            readEnvPop( conf_->optionAsString( "population_two" ), 1 );
        // translation of locations
        uint dx = conf_->optionAsInt( "grid_x" ) / 2;
        for( Population::loc_iter ii = bloc.begin();
            ii != bloc.end(); ++ii ) {
            ii->first += dx;
        }
        // set agent type to 1 and init
        for( std::vector< Agent* >::iterator ii = bux.begin(); 
            ii != bux.end(); ++ii ) {
            ( **ii ).type( 1 );
            ( **ii ).initialise();
        }
        result = new Population( xx, yy, bux, bloc,
            scalingScheme(), selectionScheme() );
    }
    return result;
}

dorsalfin::Environment*
dorsalfin::Factory::environment() {
    Environment *result = 0;
    // reading in string of bits
    std::vector< boost::dynamic_bitset<> > bux;
    readInput( bux );
    if( conf_->optionAsString( "environment" ) == "norotate" ) {
        // environment has a grid too
        uint xx = conf_->optionAsInt( "grid_x" );
        uint yy = conf_->optionAsInt( "grid_y" );
        result = new NoRotateBitEnvironment( xx, yy, 
            conf_->optionAsInt( "env_seed" ),
            conf_->optionAsDouble( "noise" ) );
        NoRotateBitEnvironment *aux =
            dynamic_cast< NoRotateBitEnvironment * >( result );
        aux->initialise( bux );
    } else if( conf_->optionAsString( "environment" ) == "noisy" ) {
        // environment has a grid too
        uint xx = conf_->optionAsInt( "grid_x" );
        uint yy = conf_->optionAsInt( "grid_y" );
        result = new NoisyBitEnvironment( xx, yy, 
            conf_->optionAsInt( "env_seed" ),
            conf_->optionAsDouble( "noise" ) );
        NoisyBitEnvironment *aux =
            dynamic_cast< NoisyBitEnvironment * >( result );
        aux->initialise( bux );
    } else if( conf_->optionAsString( "environment" ) == "subset" ) {
        // environment has a grid too
        uint xx = conf_->optionAsInt( "grid_x" );
        uint yy = conf_->optionAsInt( "grid_y" );
        result = new SubsetBitEnvironment( xx, yy, 
            conf_->optionAsInt( "env_seed" ) );
        SubsetBitEnvironment *aux =
            dynamic_cast< SubsetBitEnvironment * >( result );
        aux->initialise( bux );
    } else if( conf_->optionAsString( "environment" ) == "onebite" ) {
        // environment has a grid too
        uint xx = conf_->optionAsInt( "grid_x" );
        uint yy = conf_->optionAsInt( "grid_y" );
        result = new OneBiteBitEnvironment( xx, yy, 
            conf_->optionAsInt( "env_seed" ),
            conf_->optionAsDouble( "noise" ) );
        OneBiteBitEnvironment *aux =
            dynamic_cast< OneBiteBitEnvironment * >( result );
        aux->initialise( bux );
    } else {
        // try to read it as a file
        std::cout << "Reading environment from file" << std::endl;
        result = readEnvironment( conf_->optionAsString( "environment" ),
            bux );
    }
    // set diffusion rate
    BitEnvironment::setDiffusion( conf_->optionAsDouble( "diffusion" ) );
    return result;
}

void
dorsalfin::Factory::readInput( 
    std::vector< boost::dynamic_bitset<> > &bits ) {
    std::vector< std::string > bux( conf_->optionAsVector( "input" ) );
    // reading in string of bits
    for( std::vector< std::string >::iterator i = bux.begin(); i != bux.end();
        ++i ) {
        try {
            // decimal
            boost::dynamic_bitset<> cux( 8 * sizeof( ulong ), 
                boost::lexical_cast< ulong >( *i ) );
            bits.push_back( cux );
        } catch( boost::bad_lexical_cast &bc ) {
            // binary (and remove the '0b' header)
            boost::dynamic_bitset<> cux( i->erase( 0, 2 ) );
            bits.push_back( cux );
        }
    }
}

dorsalfin::Environment*
dorsalfin::Factory::readEnvironment( const std::string &fname, 
    std::vector< boost::dynamic_bitset<> > &bits ) {
    using namespace boost::filesystem;
    Environment *result = 0;
    std::string line, classname;
    uint xx, yy;
    ifstream infile( fname );
    // read all lines starting with a '#'
    std::getline( infile, line );
    while( line[ 0 ] == '#' ) std::getline( infile, line );
    // decipher class name, x and y
    boost::tokenizer<> tokens( line );
    boost::tokenizer<>::iterator i( tokens.begin() );
    classname = *i;
    ++i;
    xx = boost::lexical_cast< uint >( *i );
    ++i;
    yy = boost::lexical_cast< uint >( *i );
    if( classname == "NoRotateBitEnvironment" ) {
        result = new NoRotateBitEnvironment( xx, yy,
            conf_->optionAsInt( "env_seed" ),
            conf_->optionAsDouble( "noise" ) );
        dynamic_cast< BitEnvironment * >( result )->initialise( bits );
        dynamic_cast< NoRotateBitEnvironment * >( result )->readGrid( infile );
    } else if( classname == "NoisyBitEnvironment" ) {
        result = new NoisyBitEnvironment( xx, yy,
            conf_->optionAsInt( "env_seed" ),
            conf_->optionAsDouble( "noise" ) );
        dynamic_cast< BitEnvironment * >( result )->initialise( bits );
        dynamic_cast< NoisyBitEnvironment * >( result )->readGrid( infile );
    } else if( classname == "SubsetBitEnvironment" ) {
        result = new SubsetBitEnvironment( xx, yy,
            conf_->optionAsInt( "env_seed" ) );
        dynamic_cast< BitEnvironment * >( result )->initialise( bits );
        dynamic_cast< SubsetBitEnvironment * >( result )->readGrid( infile );
    } else if( classname == "OneBiteBitEnvironment" ) {
        result = new OneBiteBitEnvironment( xx, yy,
            conf_->optionAsInt( "env_seed" ),
            conf_->optionAsDouble( "noise" ) );
        dynamic_cast< BitEnvironment * >( result )->initialise( bits );
        dynamic_cast< OneBiteBitEnvironment * >( result )->readGrid( infile );
    } else {
        std::cerr << "# No environment in configuration file" << std::endl;
    }
    return result;
}

dorsalfin::ScalingScheme*
dorsalfin::Factory::scalingScheme() {
    ScalingScheme *result;
    if( conf_->optionAsString( "scaling_scheme" ) == "none" ) {
        result = new NoScaling();
        Population::threshold( THRESHOLD );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "linear" ) {
        result = new LinearScaling( conf_->optionAsDouble( "base_score" ) );
        Population::threshold( THRESHOLD );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "power" ) {
        result = new PowerScaling( conf_->optionAsDouble( "base_score" ), 
            SELECTION );
        Population::threshold( std::pow( THRESHOLD, SELECTION ) );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "exponential" ) {
        result = new ExponentialScaling( conf_->optionAsDouble( "base_score" ),
            SELECTION );
        Population::threshold( std::exp( THRESHOLD * SELECTION ) );
    } else {
        // default
        result = new NoScaling();
        Population::threshold( THRESHOLD );
    }
    return result;
}

dorsalfin::SelectionScheme*
dorsalfin::Factory::selectionScheme() {
    SelectionScheme *result;
    if( conf_->optionAsString( "selection_scheme" ) == "random" ) {
        result = new RandSelection();
    } else if( conf_->optionAsString( "selection_scheme" ) == "probalistic" ) {
        result = new ProbalisticSelection();
    } else {
        // default
        result = new RandSelection();
    }
    return result;
}


dorsalfin::GenRegAgent*
dorsalfin::Factory::readEnvAgent( const std::string &fname, int type ) {
    // Part for reading the genome
    try {
        XMLPlatformUtils::Initialize();
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Error during initialization! :\n" << message << "\n";
        XMLString::release( &message );
    }

    SAX2XMLReader *parser = XMLReaderFactory::createXMLReader();
    parser->setFeature( XMLUni::fgSAX2CoreValidation, false );
    parser->setFeature( XMLUni::fgSAX2CoreNameSpaces, false );
    // AgentReader that reads genome and network
    EnvAgentReader *doc_handler = new EnvAgentReader( this, type );
    ErrorHandler *error_handler = (ErrorHandler*) doc_handler;
    parser->setContentHandler( doc_handler );
    parser->setErrorHandler( error_handler );    

    try {
        parser->parse( fname.c_str() );
    } catch( const XMLException &toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( const SAXParseException &toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( boost::bad_lexical_cast &toCatch ) {
        std::cout << "Bad lexical cast: " << toCatch.what() << "\n";
    } catch( ... ) {
        std::cout << "Unexpected Exception \n" ;
    }

    // get that agent
    GenRegAgent *result = 0;
    if( !doc_handler->sawErrors() ) {
        result = doc_handler->agent();
    }
    doc_handler->resetDocument();
    // delete parser before calling Terminate
    delete parser;
    delete doc_handler;
    // and call Terminate
    XMLPlatformUtils::Terminate();
    
    // Now we need to read the network from a dot file...
    return result;
}

boost::tuple< std::vector< dorsalfin::Agent* >, 
    std::vector< dorsalfin::location > >
dorsalfin::Factory::readEnvPop( const std::string &fname, int type ) {
    // Read genomes
    try {
        XMLPlatformUtils::Initialize();
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Error during initialization! :\n" << message << "\n";
        XMLString::release( &message );
    }
            
    SAX2XMLReader *parser = XMLReaderFactory::createXMLReader();
    parser->setFeature( XMLUni::fgSAX2CoreValidation, false );
    parser->setFeature( XMLUni::fgSAX2CoreNameSpaces, false );

    // not giving type, coz it's read from file
    EnvPopReader *doc_handler = new EnvPopReader( this, type );
    ErrorHandler *error_handler = (ErrorHandler*) doc_handler;
    parser->setContentHandler( doc_handler );
    parser->setErrorHandler( error_handler );    

    try {
        parser->parse( fname.c_str() );
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( const SAXParseException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( boost::bad_lexical_cast &toCatch ) {
        std::cout << "Bad lexical cast: " << toCatch.what() << "\n";
    } catch( ... ) {
        std::cout << "Unexpected Exception \n" ;
    }

    // get that agent population
    std::vector< Agent* > result;
    std::vector< location > loc;
    if( !doc_handler->sawErrors() ) {
        boost::tie( result, loc ) = doc_handler->population();
    }
    doc_handler->resetDocument();

    // delete parser before calling Terminate
    delete parser;
    delete doc_handler;
    // and call Terminate
    XMLPlatformUtils::Terminate();

    // still need to read networks from dot file
    network( type );
    return boost::make_tuple( result, loc );
}

void
dorsalfin::Factory::readAgentConfigurations() {
    for( int i = 0; i != conf_->optionAsInt( "nr_agent_type" ); ++i ) {
        std::string fname =
            "agent" + boost::lexical_cast< std::string >( i ) + ".cfg";
        conf_->parseAgentFile( fname );
        setConfiguration( i );
    }
}

void
dorsalfin::Factory::setConfiguration( int type ) {
    // pre: conf_ is initialised
    SimpleAgent::birthRate( conf_->optionAsDouble( "birth_rate", type ) );
    SimpleAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );
    
    EnvAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );
    EnvAgent::maxGenomeSize( conf_->optionAsInt( "max_genome_size", type ) );
    EnvAgent::maxTposons( conf_->optionAsInt( "max_tposons", type ) );
    EnvAgent::genomePenaltyRate( 
        conf_->optionAsDouble( "genome_size_penalty", type ) );
    EnvAgent::retroposonPenaltyRate(
        conf_->optionAsDouble( "tposons_penalty", type ) );

    Chromosome::copyGeneRate( conf_->optionAsDouble( "cp_gene", type ) );
    Chromosome::removeGeneRate( conf_->optionAsDouble( "rm_gene", type ) );
    Chromosome::thresholdGeneRate( conf_->optionAsDouble( "thr_gene", type ) );
    Chromosome::tagGeneRate( conf_->optionAsDouble( "tag_gene", type ) );
    Chromosome::stateRate( conf_->optionAsDouble( "state_gene", type ) );
    
    Chromosome::copyRetroposonRate( conf_->optionAsDouble( "cp_tp", type ) );
    Chromosome::removeRetroposonRate( conf_->optionAsDouble( "rm_tp", type ) );
    Chromosome::newRetroposonRate( conf_->optionAsDouble( "new_tp", type ) );
    Chromosome::removeRepeatRate( conf_->optionAsDouble( "rm_ltr", type ) );
    Chromosome::recombinationRate( 
        conf_->optionAsDouble( "dsb_recombination", type ) );
    
    Chromosome::weightInteractionRate( 
        conf_->optionAsDouble( "weight_ia", type ) );
    Chromosome::copyInteractionRate( conf_->optionAsDouble( "cp_ia", type ) );
    Chromosome::removeInteractionRate( conf_->optionAsDouble( "rm_ia", type ) );
    Chromosome::newInteractionRate( conf_->optionAsDouble( "new_ia", type ) );
    Chromosome::tagInteractionRate( conf_->optionAsDouble( "tag_ia", type ) );
}

std::vector< std::list< dorsalfin::ChromosomeElement* > > *
dorsalfin::Factory::chromosomeParts( int type ) {
    /*// # transfacs
    int n_tf = conf_->optionAsInt( "transfactors", type );*/
    // # genes
    int n_ds = conf_->optionAsInt( "genes", type );
    // # retroposons ( inbetween other genes )
    int n_tp = conf_->optionAsInt( "tposons", type );
    // # single repeat elements ( inbetween genes etc )
    int n_rp = conf_->optionAsInt( "repeats", type );

    // create chromosome parts
    std::vector< std::list< ChromosomeElement* > > *result = 
        new std::vector< std::list< ChromosomeElement* > >();
    /*// create transcription factors
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_tf; ++i ) {
        result->back().push_back( transfactor() );
    }*/
    // create genes
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_ds; ++i ) {
        result->back().push_back( gene() );
    }
    // create retroposons
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_tp; ++i ) {
        result->back().push_back( retroposon() );
    }
    // create single repeats
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_rp; ++i ) {
        result->back().push_back( repeat() );
    }
    return result;
}

std::list< dorsalfin::ChromosomeElement* > *
dorsalfin::Factory::organiseChromosome( 
    std::vector< std::list< ChromosomeElement* > > *parts, double organise ) {
    typedef std::vector< std::list< ChromosomeElement* > >::iterator aux_iter;
    std::list< ChromosomeElement* > *result =
        new std::list< ChromosomeElement* >();
    
    // paste perfect & shuffle a bit
    // last 2 elements of 'parts' are retroposon and repeat elements
    std::list< ChromosomeElement* > ltr = parts->back();
    parts->pop_back();
    std::list< ChromosomeElement* > rp = parts->back();
    parts->pop_back();
    // add repeats
    int aux = 2 * parts->size();
    Chromosome::ce_iter i = ltr.begin();
    while( i != ltr.end() ) {
        int bux = static_cast< int >( aux * uniform() );
        if( bux % 2 == 0 ) {
            ( *parts )[ bux / 2 ].push_back( *i );
        } else {
            ( *parts )[ bux / 2 ].push_front( *i );
        }
        ++i;
    }
    // and retroposons
    i = rp.begin();
    while( i != rp.end() ) {
        int bux = static_cast< int >( aux * uniform() );
        if( bux % 2 == 0 ) {
            ( *parts )[ bux / 2 ].push_back( *i );
        } else {
            ( *parts )[ bux / 2 ].push_front( *i );
        }
        ++i;
    }
    // glue together
    for( aux_iter i = parts->begin(); i != parts->end(); ++i ) {
        result->splice( result->end(), *i );
    }
    // shuffle a little bit...
    aux = result->size();
    int nr_elem = static_cast< int >( 0.5 + ( 1.0 - organise ) * aux );
    if(  nr_elem == aux ) {
        // most annoying, copy to vector, from vector
        std::vector< ChromosomeElement* > bux( result->begin(), result->end() );
        std::random_shuffle( bux.begin(), bux.end(), rand_range< int > );
        std::copy( bux.begin(), bux.end(), result->begin() );
    } else {
        // naive implementation
        // take out elements
        std::list< ChromosomeElement* > bux;
        for( int i = 0; i < nr_elem; ++i ) {
            Chromosome::ce_iter cux = boost::next( result->begin(), 
                static_cast< int >( aux * uniform() ) );
            bux.push_back( *cux );
            result->erase( cux );
            --aux;
        }
        // put them back 
        while( !bux.empty() ) {
            result->insert( 
                boost::next( result->begin(), 
                    static_cast< int >( aux * uniform() ) ), bux.back() );
            bux.pop_back();
            ++aux;
        }
    }
    delete parts;
    return result;
}

//
// Queue of events
//
std::priority_queue< dorsalfin::Event* >*
dorsalfin::Factory::eventQueue() {
    std::priority_queue< Event* > *result = new std::priority_queue< Event* >();
    std::vector< std::string > aux( conf_->optionAsVector( "event" ) );
    // do we have any events?
    if( aux.empty() ) return result;
    // first token gives name of event
    if( aux.front() == "zerorates" ) {
        // second gives time
        result->push( 
            new ZeroRatesEvent( boost::lexical_cast< long >( aux.back() ) ) );
    } else {
        std::cerr << "Unknown event." << std::endl;
    }
    return result; 
}
