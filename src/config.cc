//
// Implementation of the configuration options.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "config.hh"

dorsalfin::Config::Config() 
    : general_( "General options" ), conf_( "Configuration" ), 
      cmdline_(), collect_( "Data logging" ), agent_( "Agents" ),
      agent_maps_() {
    // quickly intro the standard config file
    stdcfg_ = "luck.cfg";

    // 1st cmd line options, 2nd config file
    general_.add_options()
        ( "config,c", bo_po::value< std::string >()->default_value( stdcfg_ ),
          "use config file arg" )
        ( "runs,r", bo_po::value< int >()->default_value( 1 ), 
          "# simulation runs" )
        ( "overview", "print current configuration" )
        ( "version,v", "print version string" )
        ( "help,h", "produce help message" );

    conf_.add_options()
        ( "grid_x", bo_po::value< int >()->default_value( 25 ),
          "width of population grid" )
        ( "grid_y", bo_po::value< int >()->default_value( 25 ),
          "height of population grid" )
        ( "shuffle", bo_po::value< std::string >()->default_value( "false" ),
          "shuffle the grid" )
        ( "sum_fitness_threshold", 
          bo_po::value< double >()->default_value( THRESHOLD ),
          "threshold for probalistic reproduction [ 0.0, 8.0 )" )
        ( "base_score", bo_po::value< double >()->default_value( 0.2 ),
          "base score for linear scaling if all scores are equal" )
        ( "scaling_scheme", 
          bo_po::value< std::string >()->default_value( "power" ), 
          "agent score scaling ( none, linear, power, exponential )" )
        ( "selection_scheme", 
          bo_po::value< std::string >()->default_value( "probalistic" ),
          "agent selection ( random, probalistic )" )
        ( "nr_agent_type", bo_po::value< int >()->default_value( 1 ),
          "# of different kinds of agents" )
        ( "agent_placement", 
          bo_po::value< std::string >()->default_value( "random" ),
          "placement of the different agent types ( random, patch )" )
        ( "population_one", bo_po::value< std::string >(),
          "if present, read population from file" )
        ( "population_two", bo_po::value< std::string >(),
          "if present, read 2nd population from file" )
        ( "population_start", 
          bo_po::value< std::string >()->default_value( "random" ), 
          "start population with(out) diversity" )
        ( "environment", 
          bo_po::value< std::string >()->default_value( "bitset" ), 
          "environment [bitset, noisy, norotate]" )
        ( "noise", bo_po::value< double >()->default_value( 1e-3 ),
          "noise in the environment" )
        ( "diffusion", bo_po::value< double >()->default_value( 0 ),
          "diffusion freq in the environment" )
        ( "init_seed", bo_po::value< int >()->default_value( 2 ),
          "random number generator seed of initialisation" )
        ( "random_seed", bo_po::value< int >()->default_value( 2 ),
          "random number generator seed of population" )
        ( "env_seed", bo_po::value< int >()->default_value( 2 ),
          "random number generator seed of environment" )
        ( "end_time", bo_po::value< long >()->default_value( 100 ),
          "# timesteps to run" )
        ( "input", bo_po::value< std::string >()->default_value( "0" ),
          "initial input of environment" )
        ( "event", bo_po::value< std::string >()->default_value( "" ),
          "events that may occur during the simulation" );

          
    collect_.add_options()
        ( "log_path", bo_po::value< std::string >()->default_value( "~/tmp" ),
          "path where collected data is written to" )
        ( "log_mutations_csv", bo_po::value< std::string >(),
          "# dsbs, gene cp/rm" )
        ( "log_pop_grid_csv", bo_po::value< std::string >(),
          "population grid pathname" ) 
        ( "log_env_grid_csv", bo_po::value< std::string >(),
          "environment grid pathname" )
        ( "log_pop_gridline_csv", bo_po::value< std::string >(),
          "population grid line filename" )
        ( "log_env_gridline_csv", bo_po::value< std::string >(),
          "environment grid line filename" )
        ( "log_scores_csv", bo_po::value< std::string >(),
          "fitness score in csv filename" )
        ( "log_distances_csv", bo_po::value< std::string >(),
          "geno distances in csv filename" )
        ( "log_distances_histo_csv", bo_po::value< std::string >(),
          "geno distances histograms in csv filename" )
        ( "log_lengths_histo_csv", bo_po::value< std::string >(),
          "geno lengths histograms in csv filename" )
        ( "log_pop_bitset_csv", bo_po::value< std::string >(),
          "bitsets produced by agents in csv filename" )
        ( "log_env_bitset_csv", bo_po::value< std::string >(),
          "bitsets in environment in csv filename" )
        ( "log_population_csv", bo_po::value< std::string >(),
          "population sizes in csv filename" )
        ( "log_ancestors_csv", bo_po::value< std::string >(),
          "ancestor tracing in csv filename" )
        ( "log_agent_trace_xml", bo_po::value< std::string >(),
          "agent-trace-in-xml pathname" )
        ( "agent_trace_source_csv", bo_po::value< std::string >(),
          "source trace in csv" )
        ( "log_genomes_xml", bo_po::value< std::string >(),
          "genomes-in-xml pathname" )
        ( "log_sample_genomes_xml", bo_po::value< std::string >(),
          "genomes-in-xml pathname" )
        ( "log_genes_csv", bo_po::value< std::string >(),
          "gene copy numbers in csv filename" )
        ( "log_active_genes_csv", bo_po::value< std::string >(),
          "active genes in csv filename" )
        ( "log_bindingsites_csv", bo_po::value< std::string >(),
          "bindingsite copy numbers in csv filename" )
        ( "log_indegree_csv", bo_po::value< std::string >(),
          "in degree of genes in csv filename" )
        ( "log_outdegree_csv", bo_po::value< std::string >(),
          "out degree of genes in csv filename" )
        ( "log_retroposons_csv", bo_po::value< std::string >(),
          "nr of repeats and retros in csv filename" )
        ( "log_period", bo_po::value< long >()->default_value( 1 ),
          "frequency of writing to file (once every..)" )
        ( "log_period_stats", bo_po::value< long >()->default_value( 1 ),
          "frequency of writing to gene and bsite counts file (once every..)" )
        ( "log_period_xml", bo_po::value< long >()->default_value( 1 ),
          "frequency of writing to xml file (once every..)" );

    agent_.add_options()
        ( "init_nr_agents", bo_po::value< int >()->default_value( 1 ), 
          "initial # agents in population" )
        ( "cp_gene", bo_po::value< double >()->default_value( 1e-4 ),
          "rate of gene duplication" )
        ( "rm_gene", bo_po::value< double >()->default_value( 1e-4 ),
          "rate of removing an up and downstream" )
        ( "tag_gene", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of interaction targets" )
        ( "thr_gene", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of threshold per gene" )
        ( "cp_tp", bo_po::value< double >()->default_value( 1e-4 ),
          "rate of retroposon activity"  )
        ( "rm_ltr", bo_po::value< double >()->default_value( 1e-5 ),
          "rate of repeat removal" )
        ( "rm_tp", bo_po::value< double >()->default_value( 1e-4 ),
          "rate of reciprocal recombination of retroposons" )
        ( "new_tp", bo_po::value< double >()->default_value( 1e-6 ),
          "rate of new retroposons influx" )
        ( "dsb_recombination", bo_po::value< double >()->default_value( 1e-2 ),
          "rate of DSB occurrence * rate of faulty NHEJ" )
        ( "tposons", bo_po::value< int >()->default_value( 10 ),
          "# retrotransposons per chromosome" )
        ( "repeats", bo_po::value< int >()->default_value( 0 ),
          "# single repeats per chromosome" )
        ( "genes", bo_po::value< int >()->default_value( 20 ),
          "# genes per chromosome" )
        ( "chromos", bo_po::value< int >()->default_value( 1 ),
          "# chromosomes per genome" )
        ( "organised", bo_po::value< double >()->default_value( 0.0 ),
          "how perfect is a chromosome" )
        ( "birth_rate", bo_po::value< double >()->default_value( 0.2 ),
          "birth rate of simple agent" )
        ( "death_rate", bo_po::value< double >()->default_value( 0.1 ),
          "death rate of agents" )
        ( "agent", bo_po::value< std::string >(),
          "type of agent ( simple, net, delta )" )
        ( "genome_size_penalty", bo_po::value< double >()->default_value( 1.0 ),
          "penalty for genome size conservation" )
        ( "max_genome_size", bo_po::value< int >()->default_value( 400 ),
          "maximum size of genome" )
        ( "tposons_penalty", bo_po::value< double >()->default_value( 1.0 ),
          "penalty for having too many retroposons" )
        ( "max_tposons", bo_po::value< int >()->default_value( 50 ),
          "maximum # retrotransposons allowed" )
        ( "network", bo_po::value< std::string >()->default_value( "net.dot" ),
          "reference network" )
        ( "init_propagate", bo_po::value< int >()->default_value( 2 ),
          "init time to propagate network, before output" )
        ( "max_propagate", bo_po::value< int >()->default_value( 10 ),
          "max time to calculate stable state of network" )
        ( "seq_propagate", 
          bo_po::value< std::string >()->default_value( "false" ),
          "sequential or parallel updating of the network" )
        ( "weight_ia", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of interaction weights" )
        ( "cp_ia", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of copying interactions" )
        ( "rm_ia", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of removing interactions" )
        ( "new_ia", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of neogenesis" )
        ( "tag_ia", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation rate of interaction targets" )
        ( "state_gene", bo_po::value< double >()->default_value( 1e-4 ),
          "mutation or noise rate of state switching" );

/*        ( "max_distance", bo_po::value< int >()->default_value( 4 ),
          "maximum distance from ideal genotype having a nonzero fitness" )
*/

    cmdline_.add( general_ ).add( conf_ ).add( collect_ );
    // dynamically create variables map, coz we want to create a new one from
    // time to time
    var_map_ = new bo_po::variables_map();
}

dorsalfin::Config::~Config() {
    if( var_map_ != 0 ) {
        delete var_map_;
    }
    // and delete agent configs
    agent_maps_.clear();
}

void
dorsalfin::Config::reset() {
    if( var_map_ != 0 ) {
        delete var_map_;
    }
    var_map_ = new bo_po::variables_map();
    //agent_maps_.clear();
}

void
dorsalfin::Config::parseCmdLine( int argc, char **argv ) {
    bo_po::store( 
            bo_po::parse_command_line( argc, argv, cmdline_ ), *var_map_ );
    bo_po::notify( *var_map_ );
}

void 
dorsalfin::Config::parseFile() {
    std::string fname = ( *var_map_ )[ "config" ].as< std::string >(); 
    doParseFile( fname , cmdline_, *var_map_ );
}

void 
dorsalfin::Config::parseFile( std::string fname ) {
    doParseFile( fname , cmdline_, *var_map_ );
}

void 
dorsalfin::Config::parseAgentFile( std::string fname ) {
    agent_maps_.push_back( bo_po::variables_map() );
    doParseFile( fname, agent_, agent_maps_.back() );
}

int 
dorsalfin::Config::optionAsInt( const std::string &s ) {
    return ( *var_map_ )[ s ].as< int >();
}

int 
dorsalfin::Config::optionAsInt( const std::string &s, int agent ) const {
    return agent_maps_[ agent ][ s ].as< int >();
}

long 
dorsalfin::Config::optionAsLong( const std::string &s ) {
    return ( *var_map_ )[ s ].as< long >();
}

long 
dorsalfin::Config::optionAsLong( const std::string &s, int agent ) const {
    return agent_maps_[ agent ][ s ].as< long >();
}

double 
dorsalfin::Config::optionAsDouble( const std::string &s ) {
    return ( *var_map_ )[ s ].as< double >();
}

double 
dorsalfin::Config::optionAsDouble( const std::string &s, int agent ) const {
    return agent_maps_[ agent ][ s ].as< double >();
}

std::string 
dorsalfin::Config::optionAsString( const std::string &s ) {
    return ( *var_map_ )[ s ].as< std::string >();
}

std::string 
dorsalfin::Config::optionAsString( const std::string &s, int agent ) const {
    return agent_maps_[ agent ][ s ].as< std::string >();
}

std::vector< std::string >
dorsalfin::Config::optionAsVector( const std::string &s ) {
    std::vector< std::string > result;    
    boost::tokenizer<> tok( ( *var_map_ )[ s ].as< std::string >() );
    for( boost::tokenizer<>::iterator i = tok.begin(); i != tok.end(); ++i ) {
        result.push_back( *i );
    }
    return result;
}

void
dorsalfin::Config::help( std::ostream &os ) const {
    os << cmdline_ << "\n";
}

void
dorsalfin::Config::version( std::ostream &os ) const {
    os << "Dorsalfin version " << dorsalfin::VERSION << "\n";
}

void
dorsalfin::Config::overview( std::ostream &os ) const {
    // First print version
    version( os );
    // slightly dirty programming..
    int aux = os.width();
    int bux = 20;
    for( var_iter i = var_map_->begin(); i != var_map_->end(); ++i ) {
        os.width( bux );
        os << std::left << i->first;
        os.width( aux );
        os << ": ";
        
        try {
            os << ( i->second ).as< int >();
        } catch( boost::bad_any_cast ) {}
        try {
            os << ( i->second ).as< long >();
        } catch( boost::bad_any_cast ) {}
        try {
            os << ( i->second ).as< std::string >();
        } catch( boost::bad_any_cast ) {}
        try {
            os << ( i->second ).as< double >();
        } catch( boost::bad_any_cast ) {}
        os << "\n";
    }
    os.flush();
}

void
dorsalfin::Config::write( std::ostream &os ) const {
    // quite similar to overview
    os << "# simulation parameters\n";
    os << "version = " << dorsalfin::VERSION << "\n";
    for( var_iter i = var_map_->begin(); i != var_map_->end(); ++i ) {
        os << i->first << " = ";        
        try {
            os << ( i->second ).as< int >();
        } catch( boost::bad_any_cast ) {}
        try {
            os << ( i->second ).as< long >();
        } catch( boost::bad_any_cast ) {}
        try {
            os << ( i->second ).as< std::string >();
        } catch( boost::bad_any_cast ) {}
        try {
            os << ( i->second ).as< double >();
        } catch( boost::bad_any_cast ) {}
        os << "\n";
    }
    // one extra
    os << "selection = " << SELECTION << "\nthreshold = " << THRESHOLD << "\n";
    
    // and now the agent configs
    os << "\n# agent configuration\n";
    typedef std::vector< bo_po::variables_map >::const_iterator ag_iter;
    for( ag_iter j = agent_maps_.begin(); j != agent_maps_.end(); ++j ) {
        for( var_iter i = j->begin(); i != j->end(); ++i ) {
            os << i->first << " = ";        
            try {
                os << ( i->second ).as< int >();
            } catch( boost::bad_any_cast ) {}
            try {
                os << ( i->second ).as< long >();
            } catch( boost::bad_any_cast ) {}
            try {
                os << ( i->second ).as< std::string >();
            } catch( boost::bad_any_cast ) {}
            try {
                os << ( i->second ).as< double >();
            } catch( boost::bad_any_cast ) {}
            os << "\n";
        }
        os << "\n";
    }
    os.flush();
}

void
dorsalfin::Config::doParseFile( std::string &fname, 
    bo_po::options_description &options, bo_po::variables_map &vmap ) {
    boost::filesystem::fstream file;
    file.open( fname, boost::filesystem::fstream::in );

    if( file.is_open() ) {
        bo_po::store( bo_po::parse_config_file( file, options ), vmap );
        bo_po::notify( vmap );
    } else {
        std::string msg = "could not open any config file named ";
        msg += fname;
        throw msg.c_str();
    }
    file.close();
}
