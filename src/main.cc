//
// Starting point of the fluke program.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "defs.hh"
#include "fluke.hh"

#ifdef DEBUG
#include <boost/timer.hpp>
#endif

using namespace dorsalfin;

// globals...
base_generator_type dorsalfin::generator( 18 );
uniform_gen_type dorsalfin::uniform( dorsalfin::generator );

int
main( int argc, char **argv ) {
    int result = 1;
    try {
#ifdef DEBUG
        boost::timer tt;
#endif
        Fluke f( argc, argv );
        result = f.run();
#ifdef DEBUG
        std::cout << "Timing: " << tt.elapsed() << " seconds\n";
#endif
    } catch( char *e ) {
        std::cout << "Exception: " << e << std::endl;
    } catch( char const *e ) {
        std::cout << "Exception: " << e << std::endl;
    } catch( std::exception &e ) {
        std::cout << "Exception: " << e.what() << std::endl;
    }
    return result;
}

