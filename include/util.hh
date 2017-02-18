//
// Some extra functions missing in STL and Boost.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_UTIL_
#define _DORSALFIN_UTIL_

namespace dorsalfin {

    /// Reverse (or swap) the two items of a pair
    template< typename T1, typename T2 > std::pair< T2, T1 >
    reverse( std::pair< T1, T2 > p ) {
        return std::make_pair( p.second, p.first );
    }
    
    /// A variant of the STL function \c std::for_each. 
    /// \pre \c first2 is not in the range [\c first, \c last)
    template< typename In, typename In2, typename BinOp > BinOp
    for_each( In first, In last, In2 first2, BinOp op ) {
        while( first != last ) op( *first++, *first2++ );
        return op;
    }
    
    /// Copy if a predicate holds. A function clearly missing from the STL.
    /// \pre \c res is not in the range [\c first, \c last)
    template< typename In, typename Out, typename Pred > Out
    copy_if( In first, In last, Out res, Pred p ) {
        while( first != last ) {
            if( p( *first ) ) *res++ = *first;
            ++first;
        }
        return res;
    }

    /// Indirect erase. Normally the STL function family \c erase deletes 
    /// only the object in the container. In \b fluke most of the containers 
    /// hold pointers to objects. \c smart_erase erases the object pointed to
    /// as well.
    template< typename Container, typename For > For
    smart_erase( Container &x, For first ) {
        delete *first;
        return x.erase( first );
    }

    /// Indirect range erase. Both the range and the objects pointed to are
    /// erased.
    /// See also smart_erase( Container, For )
    template< typename Container, typename For > For
    smart_erase( Container &x, For first, For last ) {
        For i = first;
        while( i != last ) {
            delete *i;
            ++i;
        }
        return x.erase( first, last );
    }
    
    /// Numerically safe floating-point comparison.
    template< typename T > bool
    close_to( T f, T g ) {
        return std::abs( f - g ) < std::numeric_limits< T >::epsilon();
    }

    /// Mean
    template< typename In > double
    mean( In first, In last ) {
        double mm = 0.0;
        mm = std::accumulate( first, last, mm );
        return mm / std::distance( first, last );
    }

    /// Variance
    template< typename In > double
    variance( In first, In last, double mean ) {
        double aux = std::distance( first, last );
        double bux = 0.0;
        while( first != last ) {
            bux += ( *first - mean ) * ( *first - mean );
            ++first;
        }
        return bux / ( aux - 1 );
    }
    
    /// Absolute
    template< typename T >
    struct absolute : public std::unary_function< T, T > {
        T operator()( const T &x ) const 
        { return x < 0? -x: x; }
    };

    /// Retrieve random element
    template< typename For, typename RandomGenerator > For
    random_element( For first, For last, RandomGenerator &rangen ) {
        int aux = std::distance( first, last );
        return boost::next( first, rangen( aux ) );
    }
    
    /// Swap two bits
    /// pre: i, j < b.size()
    template< typename Block, typename Allocator > void
    swap_bits( boost::dynamic_bitset< Block, Allocator > &b, uint i, uint j ) {
        bool tmp = b[ i ];
        b[ i ] = b[ j ];
        b[ j ] = tmp;
    }
    
    /// Rotate bitset to the right
    /// pre: offset < b.size()
    ///
    /// Naive implementation, counting on loop unrolling etc to speed this up
    template< typename Block, typename Allocator > void
    rotate( boost::dynamic_bitset< Block, Allocator > &b, uint offset ) {
        if( offset == 0 or offset == b.size() ) return;
        boost::dynamic_bitset< Block, Allocator > a( b );
        a <<= ( b.size() - offset );
        b = a | ( b >> offset );
    }
    
    /// Reverse boost::dynamic_bitset
    template< typename Block, typename Allocator > void
    reverse_bitset( boost::dynamic_bitset< Block, Allocator > &b ) {
        uint mid = b.size() / 2;
        uint tmp = b.size() - 1;
        for( uint i = 0; i != mid; ++i ) {
            swap_bits( b, i, tmp - i );
        }
    }
}
#endif

