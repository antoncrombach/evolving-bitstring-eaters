//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_RETROPOSON_H_
#define _DORSALFIN_RETROPOSON_H_

#include "chromelement.hh"

namespace dorsalfin {

    /// \class Retroposon
    /// \brief Retrotransposable element.
    ///
    /// Retroposons are a known mechanism for gene duplication and for the 
    /// occurrence of double-strand breaks (DSBs). 
    /// A retroRetroposon is always flanked by long terminal repeats, 
    /// represented by elements of \c Repeat.
    class Retroposon : public ChromosomeElement {
        public:
        /// Constructor.
        Retroposon() : ChromosomeElement() {}
        /// Constructor. The \c int is a unique identifier.
        Retroposon( uint t ) : tag_( t ) {}
        /// Copy constructor.
        Retroposon( const Retroposon &rp )
        { Retroposon::copy( rp ); }
        /// Destructor.
        virtual ~Retroposon() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another retroposon into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return to pool.
        virtual bool toPool();
        
        /// Write a html representation.
        virtual void writeXml( std::ostream &os ) const;
        /// Write a string representation
        virtual std::string asString() const;
        /// Write a label for graphs.
        virtual std::string asLabel() const;

        private:
        uint tag_;
    };

    inline std::string Retroposon::asLabel() const
    { return std::string( "R" ); }
    
    inline std::string Retroposon::asString() const
    { return boost::lexical_cast< std::string >( tag_ ); }
}
#endif

