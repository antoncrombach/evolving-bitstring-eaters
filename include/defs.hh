//
// Definitions of globals needed by entire project.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_DEFS_H_
#define _DORSALFIN_DEFS_H_

// We do not want multi_array range checking in production code
#ifndef DEBUG
#define BOOST_DISABLE_ASSERTS
#endif

#include <unistd.h>

#include <map>
#include <list>
#include <cmath>
#include <queue>
#include <string>
#include <stack>
#include <vector>
#include <cassert>
#include <limits>
#include <utility>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <ext/algorithm>
#include <ext/functional>
#include <ext/numeric>

#include <boost/regex.hpp>
#include <boost/utility.hpp>
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/lagged_fibonacci.hpp>

#include <boost/graph/random.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/vector_property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "util.hh"
#include "pool.hh"
#include "xml_inout.hh"

//
// Memory policy for arguments and return values: 
// - expect aliases (pointers and references) as arguments. Methods need to 
//   copy things themselves.
// - expect aliases as return values.
//
// Option: do not return objects, instead pass them as a reference argument.
// Removes the hazard of implicitly copying objects.

#ifdef DEBUG
using std::cout;
using std::endl;
#endif

// Random number generators
typedef boost::lagged_fibonacci607 base_generator_type;
typedef boost::variate_generator< base_generator_type, boost::uniform_int<> > 
    randrange_gen_type;
typedef boost::uniform_01< base_generator_type > uniform_gen_type;

/// \namespace dorsalfin The project is placed in the namespace \c dorsalfin.
/// Another part of a dolphin. Other namespaces used (only explicitly though)
/// are \c std and \c boost.
///
/// Some design decisions that have been made (and may be reconsidered) are 
/// listed below.
/// - the id of an agent is (time of birth, x, y, i), where (x, y) is its grid
///   location and i is the index for multiple divisions in 1 timestep
/// - distinction between lifetime and replication mutations is not in place
/// - using xerces-c SAX2 parser for xml access of XML agent and population
///   reader
/// - xml is only nice for saving a simulation or whole genomes 
/// - extra options needed for g++ to be able to write >2Gb files
///   (-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE64_SOURCE)    
///   (fixed in libstd6.0+)
///
/// Compiler directives are:
/// - SENSOR; to get an environmental sensor going,
/// - PERTURB; to have code for perturbations of state during the lifetime,
/// - SEQUENTIAL; to be able to choose between parallel or sequential updating.
///
namespace dorsalfin {

    /// Version number of the application in the format "a.b.c". \c a is 
    /// increased when either major changes are made to the program or a fully
    /// stable, 'public' release is made. \c b signals the inclusion of new
    /// features in the current release. \c c is used for series of bugfixes.
    const std::string VERSION = "dorsalfin 1.9.0";

    /// location
    typedef std::pair< int, int > location;

    /// constants
    const double PI = 3.14159265;
    // Selection pressure
    const double SELECTION = 1.0;
    // The threshold given here still needs to be rescaled according to the 
    // selection coefficient
    // Note: the actual THRESHOLD is computed in the factory
    const double THRESHOLD = 1.0;

    // all the classes declared
    class ChromosomeElement;
    class Gene;
    /*class TransFactor;*/
    class Interaction;
    class Retroposon;
    class Repeat;
    class Chromosome;
    class Genome;

    class RefNode;
    class Network;
    
    class Agent;
    class AgentTag;
    class DuoAgent;
    class SimpleAgent;
    class GenRegAgent;
    class EnvAgent;

    class Environment;
    class BitEnvironment;
    class NoisyBitEnvironment;
    class NoRotateBitEnvironment;
    
    class Population;
    class ScalingScheme;
    class NoScaling;
    class LinearScaling;
    class SelectionScheme;
    class RandSelection;
    class ProbalisticSelection;

    // rethink with concept of owner and being observed uncoupled
    class Observer;
    class LogObserver;
    class AsyncLogObserver;
    class Subject;

    class LogCsvMutations;
    class LogCsvScores;
    class LogCsvDistances;
    class LogCsvHistoDistances;
    class LogCsvHistoLengths;
    
    class LogCsvGenes;
    class LogCsvBindingSites;
    class LogCsvInDegree;
    class LogCsvOutDegree;
    class LogCsvRetroposons;
    
    class LogCsvAncestors;
    class LogXmlAgentTrace;
    class LogXmlGenomes;
    class LogXmlSampleGenomes;
    
    class LogPopulationSize;
    class LogCsvPopulationBitset;
    class LogCsvPopulationGrid;
    class LogCsvEnvironmentBitset;
    class LogCsvEnvironmentGrid;
    
    class Factory;
    class ObserverManager;
    class StreamManager;
    class Model;
    class Config;
    class Fluke;

    class Event;
    class ZeroRatesEvent;

    /// Random number generator. It is a global object as it is used 
    /// throughout the entire program. It needs to be used in conjunction
    /// with a distribution.
    extern base_generator_type generator;
    /// Uniform random numbers [0,1). It is a global object to provide easy
    /// access in the entire program.
    extern uniform_gen_type uniform;
    
    /// Generate random number of type T from the interval \f$[0,n)\f$.
    template< typename T > T
    rand_range( T n ) {
        return static_cast< T >( n * uniform() );
    }
}
#endif

