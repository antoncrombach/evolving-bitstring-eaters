//
// Xml input/output base class
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_XML_INOUT_H_
#define _DORSALFIN_XML_INOUT_H_

namespace dorsalfin {
    
    /// \class XmlWriter
    /// \brief Interface for writing to/reading from an (xml) stream
    class XmlWriter {
        public:
        /// Virtual dtor.
        virtual ~XmlWriter() {};
        
        /// Signature of xml writing.
        virtual void writeXml( std::ostream & ) const = 0;
        /// Signature of xml reading.
        //virtual void readXml( std::istream & ) = 0;
        
        protected:
        /// Hidden ctor.
        XmlWriter() {};
        /// Hidden copy ctor.
        explicit XmlWriter( const XmlWriter &x ) {};
    };

/// Overloaded \c << operator for easy writing to streams.
inline std::ostream& operator<<( std::ostream& os, const XmlWriter& xio )
{ xio.writeXml( os ); return os; }

}
#endif
