//
// Configuration options of fluke.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_CONFIG_H_
#define _DORSALFIN_CONFIG_H_

#include "defs.hh"

namespace dorsalfin {

    namespace bo_po = boost::program_options;
    
    /// \class Config
    /// \brief Store the configuration settings of the simulation.
    ///
    /// The class \c Config is responsible for parsing the command line, 
    /// reading configuration files and making all the configuration settings
    /// available to the rest of the program.
    ///
    /// Different categories of configuration settings are distinguished. 
    /// General options include version number and displaying the help message.
    /// Configuration options are usually given in a file, data collection
    /// options as well. All can be supplied via the command line or in a 
    /// configuration file. The command line overrides the options given in the
    /// file.
    ///
    /// It can also provide an overview of the current settings, or a help
    /// message explaining the use of all configuration parameters.
    class Config {
        public:
            /// Configuration items iterator
            typedef bo_po::variables_map::const_iterator var_iter;
            
        public:
            /// Constructor
            Config();
            /// Destructor
            ~Config();

            /// Reset the current configuration values. Note that the 
            /// information on how many simulation runs to do is erased too
            void reset();
            /// Parse the command line for options
            void parseCmdLine( int, char** );
            /// Parse the main config file
            void parseFile();
            /// Parse given config file
            void parseFile( std::string );
            /// Parse an agent config file
            void parseAgentFile( std::string );
            
            /// Return the option as an integer
            int optionAsInt( const std::string & );
            /// Return the option as a long integer
            long optionAsLong( const std::string & );
            /// Return the option as a floating point
            double optionAsDouble( const std::string & );
            /// Return the option as a string
            std::string optionAsString( const std::string & );
            /// Return the option as a vector of integers
            std::vector< std::string > optionAsVector( const std::string & );

            /// Return the option as an integer
            int optionAsInt( const std::string &, int ) const;
            /// Return the option as a long integer
            long optionAsLong( const std::string &, int ) const;
            /// Return the option as a floating point
            double optionAsDouble( const std::string &, int ) const;
            /// Return the option as a string
            std::string optionAsString( const std::string &, int ) const;

            /// Check if option is supplied either via command line or
            /// configuration file
            bool hasOption( const std::string &s ) const;
            /// Check if the help message is needed
            bool needsHelp() const;
            /// Check if the version of fluke is needed
            bool needsVersion() const;
            /// Check if an overview of the current configuration is needed
            bool needsOverview() const;
            /// Output the help message to an \c ostream
            void help( std::ostream & ) const;
            /// Output the current version number to an \c ostream
            void version( std::ostream & ) const;
            /// Output the overview of the current configuration to 
            /// an \c ostream
            void overview( std::ostream & ) const;

            /// Write a short output of the configuration to an \c ostream
            void write( std::ostream & ) const;

        protected:
            void doParseFile( std::string &, bo_po::options_description &, 
                bo_po::variables_map & );
            
        private:
            bo_po::options_description general_;
            bo_po::options_description conf_;
            bo_po::options_description cmdline_;
            bo_po::options_description collect_;
            bo_po::options_description agent_;
            
            bo_po::variables_map *var_map_;
            std::vector< bo_po::variables_map > agent_maps_;
            std::string stdcfg_;
    };
    
    /// Overloaded \c ostream operator for outputting the configuration
    inline std::ostream& operator<<( std::ostream &os, const Config &cfg )
    { cfg.write( os ); return os; }
    
    inline bool Config::hasOption( const std::string &s ) const
    { return var_map_->count( s ) > 0; }

    inline bool Config::needsHelp() const
    { return var_map_->count( "help" ) > 0; }

    inline bool Config::needsVersion() const 
    { return var_map_->count( "version" ) > 0; }

    inline bool Config::needsOverview() const 
    { return var_map_->count( "overview" ) > 0; }
}
#endif

