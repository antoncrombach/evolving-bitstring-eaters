//
// The heart of the program 'fluke'. A very simple heart :)
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "fluke.hh"
#include "config.hh"
#include "model.hh"
#include "stream_manager.hh"


dorsalfin::Fluke::Fluke( int argc, char **argv ) 
    : model_( 0 ), config_fname_( "" ) {
    // init configuration
    config_ = new Config();
    std::cout << "Parsing command line.." << std::endl;
    config_->parseCmdLine( argc, argv );
    if( !config_->needsVersion() and !config_->needsHelp() ) {
        std::cout << "Parsing configuration file.." << std::endl;
        config_->parseFile();
        config_fname_ = config_->optionAsString( "config" );
                
        // bug: somehow the default values of the configuration are not 
        // available before a file has been parsed?

        // init globals
        // note: no calls to the random number generator may have been made at
        // this point
        generator.seed( config_->optionAsInt( "init_seed" ) );
        uniform_gen_type bux( generator );
        uniform = bux;
        
        std::string aux = config_->optionAsString( "log_path" );
        StreamManager::instance( aux );
        model_ = new Model( this );
    }
}

dorsalfin::Fluke::~Fluke() {
    if( model_ != 0 ) delete model_;
    delete config_;
}

int
dorsalfin::Fluke::run() {
    // run the program
    if( config_->needsHelp() ) {
        config_->help( std::cout );
    } else if( config_->needsVersion() ) {
        config_->version( std::cout );
    } else if( config_->needsOverview() ) {
        config_->overview( std::cout );
    } else {
        simulate();
    }
    return 0;
}

void
dorsalfin::Fluke::simulate() {
    // how many times?
    int runs = config_->optionAsInt( "runs" );
    
    // init model
    std::cout << "Building.." << std::endl;
    model_->build();
    
    // and start the simulation random seed generator
    generator.seed( config_->optionAsInt( "random_seed" ) );
    uniform_gen_type bux( generator );
    uniform = bux;

    // one simulation is always run 
    doRun();
    for( int rr = 1; rr != runs; ++rr ) {
        reconfigure( rr );
        doRun();
    }
}

void
dorsalfin::Fluke::doRun() {
    std::cout << "Initializing.." << std::endl;
    model_->initialise();
    
    StreamManager *stream = StreamManager::instance();
    // write simulation parameters to file just be4 we start running
    boost::filesystem::ofstream *aux = 
        stream->openOutFileStream( std::string( "simulation.cfg" ),
            std::fstream::out );
    *aux << *config_;
    stream->closeOutFileStream( aux );
    
    // run simulation
    std::cout << "Running.." << std::endl;
    while( not model_->hasEnded() ) {
        model_->step();
    }
    model_->finish();
}

void
dorsalfin::Fluke::reconfigure( int r ) {
    // sort-of rerunning a modified version of the constructor
    // round up
    config_->reset();

    // assuming the file is present...
    std::string bux( config_fname_ );
    std::string::size_type aux = bux.find_last_of( "." );
    bux.replace( aux - 1, 1, boost::lexical_cast< std::string >( r ) );
    config_->parseFile( bux );
    
#ifdef DEBUG
    cout << "! rebuild" << endl;
#endif
    // new output path
    StreamManager::instance()->createSimulationPath();
    // start with correct population
    model_->rebuild();
    // and start the simulation random seed generator
    generator.seed( config_->optionAsInt( "random_seed" ) );
    uniform_gen_type cux( generator );
    uniform = cux;
}

