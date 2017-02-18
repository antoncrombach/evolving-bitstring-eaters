//
// Manage file streams.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_STREAM_H_
#define _DORSALFIN_STREAM_H_

#include "defs.hh"
// 64 bits file support, for writing >2Gb files
#include <ext/stdio_filebuf.h>
#include <cstdio>


namespace dorsalfin {

    /// \class StreamManager
    /// \brief Handles the file streams of the program.
    ///
    /// In order to hide the details of communicating with the underlying
    /// filesystem and to centralise all such communication a bit more, 
    /// \c StreamManager provides a small set of methods for opening, closing
    /// filestreams and creating new paths within the data path.
    ///
    /// For the logging of data, we use a systematic way of naming. Within a
    /// user-specified directory (default \c $HOME/data) a simulation 
    /// directory is created. It is of the format \c xxx-xxx, where
    /// \f$ x \in (0,1,..,9)\f$. Each new simulation has a new number
    /// (incrementing the maximum with one) and a symbolic link, \c latest,
    /// is provided to the latest simulation.
    class StreamManager {
        public:
        /// Splitting a single string gives a vector of strings.
        typedef std::vector< std::string > split_vector_type;
        /// Iterator for a splitt string.
        typedef std::vector< std::string >::iterator split_iter;
        /// Iterators for a list of open file streams.
        typedef std::list< boost::filesystem::ifstream* >::iterator infile_iter;
        typedef std::list< boost::filesystem::ofstream* >::iterator 
            outfile_iter;

        public:
        /// Destructor.
        ~StreamManager();

        /// Open a file stream for input. It mimicks the behaviour of 
        /// the standard \c std::fstream.
        boost::filesystem::ifstream* openInFileStream( std::string,
               std::ios_base::openmode );
        /// Open a file stream for output. It mimicks the behaviour of 
        /// the standard \c std::fstream.
        boost::filesystem::ofstream* openOutFileStream( std::string,
               std::ios_base::openmode );
        /// Close a given infile stream.
        void closeInFileStream( boost::filesystem::ifstream * );
        /// Close a given file stream.
        void closeOutFileStream( boost::filesystem::ofstream * );
        /// Close all open file streams.
        void closeAllFileStreams();

        /// Open (create) a directory within the simulation directory.
        void openPath( std::string );
        /// Explicitly create a new simulation directory
        void createSimulationPath();

        public:
        /// Get that one instance that's alive
        static StreamManager * instance();
        /// Get and create one instance
        static StreamManager * instance( const std::string & );
        
        protected:
        /// Constructor.
        StreamManager();
            
        private:
        void simulationPath();
        boost::filesystem::path formatSimFolder( int );
        bool isDataFolder( boost::filesystem::directory_iterator ) const;
            
        private:
        boost::filesystem::path basepath_;
        boost::filesystem::path simulation_folder_;
        std::list< boost::filesystem::ifstream* > in_files_;
        std::list< boost::filesystem::ofstream* > out_files_;
        std::map< boost::filesystem::ofstream*, FILE* > open_files_;
        
        private:
        static StreamManager *instance_;
    };
}
#endif

