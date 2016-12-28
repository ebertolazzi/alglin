/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef BORDERED_CR_HH
#define BORDERED_CR_HH

#include "Alglin.hh"
#include "Alglin++.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace alglin {

  /*\
   |   ___             _                _    ___ ___
   |  | _ ) ___ _ _ __| |___ _ _ ___ __| |  / __| _ \
   |  | _ \/ _ \ '_/ _` / -_) '_/ -_) _` | | (__|   /
   |  |___/\___/_| \__,_\___|_| \___\__,_|  \___|_|_\
  \*/

  //! Cyclic reduction of a block bidiagonal matrix
  /*!
   * 
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   *
   */
  template <typename t_Value>
  class BorderedCR {
  private:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    BorderedCR(BorderedCR const &) ;
    BorderedCR const & operator = (BorderedCR const &) ;

  protected:

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseInteger ;

    integer nblock ; //!< total number of blocks
    integer n      ; //!< size of square blocks
    integer q      ; //!< extra BC
    integer nb     ; //!< border size

    // some derived constants
    integer nx2 ;
    integer nxn ;
    integer nxnb ;

    void
    buildT( valueConstPointer TOP,
            valueConstPointer BOTTOM,
            valuePointer      T,
            integer *         iperm ) const ;

    void
    applyT( valueConstPointer T,
            integer const *   iperm,
            valuePointer      TOP,
            valuePointer      BOTTOM,
            integer           ncol ) const ;

    void
    applyT( valueConstPointer T,
            integer const *   iperm,
            valuePointer      TOP,
            valuePointer      BOTTOM ) const ;

    Factorization<t_Value> * la_factorization ;
    Factorization<t_Value> * bb_factorization ;

    /*
    //
    //  Matrix structure
    //
    //                n * (nblock+1)
    //    ___________________^____________________
    //   /                                        \
    //     n     n     n                        n    nb  q
    //  +-----+-----+-----+----.........-----+-----+---+--+   -+
    //  |  D  |  E  |     |                  |     | B |  | n  |
    //  +-----+-----+-----+             -----+-----+---+--+    |
    //  |     |  D  |  E  |                  |     | B |  | n  |
    //  +-----+-----+-----+-----+       -----+-----+---+--+    |
    //  |     |     |  D  |  E  |            |     | B |  | n  |
    //  +-----+-----+-----+-----+       -----+-----+---+--+    |
    //  :                                                 :    |
    //  :                                                 :    |
    //  :                                                 :     > n * nblock
    //  :                                                 :    |
    //  :                                                 :    |
    //  :                        +-----+-----+-----+---+--+    |
    //  :                        |  D  |  E  |     | B |  | n  |
    //  :                        +-----+-----+-----+---+--+    |
    //  :                              |  D  |  E  | B |  | n  |
    //  +-----+-----+---......---+-----+-----+-----+---+--+   -+
    //  |  C  |  C  |            |  C  |  C  |  C  | F |  |    | nb
    //  +=====+=====+===......===+=====+=====+=====+===+==+   -+
    //  |     |                              |     |   |  |    |
    //  |  H0 |                              |  HN |Hp |Hq|    | n+q
    //  |     |                              |     |   |  |    |
    //  +-----+-----+---......---+-----+-----+-----+---+--+   -+
    //                                                   q
    */

    valuePointer H0Npq  ;
    valuePointer Bmat, Cmat, Dmat, Emat, Fmat, Tmat, Hmat, Work ;
    integer      *Tperm, *Hperm ;

  private:

    LU<t_Value>  la_lu,  bb_lu  ;
    QR<t_Value>  la_qr,  bb_qr  ;
    QRP<t_Value> la_qrp, bb_qrp ;
    SVD<t_Value> la_svd, bb_svd ;

  public:

    explicit
    BorderedCR()
    : baseValue("BorderedCR_values")
    , baseInteger("BorderedCR_integers")
    , nblock(0)
    , n(0)
    , q(0)
    , nb(0)
    , nx2(0)
    , nxn(0)
    , nxnb(0)
    , la_factorization(&la_lu)
    , bb_factorization(&bb_lu)
    , H0Npq(nullptr)
    , Bmat(nullptr)
    , Cmat(nullptr)
    , Dmat(nullptr)
    , Emat(nullptr)
    , Fmat(nullptr)
    , Tmat(nullptr)
    {}

    virtual ~BorderedCR()
    {}

    //! load matrix in the class
    void
    allocate( integer _nblock,
              integer _n,
              integer _q,
              integer _nb ) {

      nblock = _nblock ;
      n      = _n ;
      q      = _q ;
      nb     = _nb ;
      nx2    = n*2 ;
      nxn    = n*n ;
      nxnb   = n*nb ;

      integer N    = nx2+nb+q ;
      integer wnnz = max(N,2*n*max(n,nb)) ;
      integer nnz  = wnnz+N*N+nb*nb+nxnb+nblock*(2*nxnb+4*nxn)+(n+q)*(nx2+nb+q) ;
      baseValue.allocate(size_t(nnz)) ;
      baseInteger.allocate(size_t(nblock*n+N)) ;

      Bmat  = baseValue(size_t(nblock*nxnb)) ;
      Cmat  = baseValue(size_t((nblock+1)*nxnb)) ;
      Dmat  = baseValue(size_t(nblock*nxn)) ;
      Emat  = baseValue(size_t(nblock*nxn)) ;
      Fmat  = baseValue(size_t(nb*nb)) ;
      Tmat  = baseValue(size_t(2*nblock*nxn)) ;
      Hmat  = baseValue(size_t(N*N)) ;
      H0Npq = baseValue(size_t((n+q)*(nx2+nb+q))) ;

      Work  = baseValue(size_t(wnnz)) ;

      Tperm = baseInteger(size_t(nblock*n)) ;
      Hperm = baseInteger(size_t(N)) ;
    }

    // filling bidiagonal part of the matrix
    void
    loadD( integer nbl, valueConstPointer D, integer ldD )
    { gecopy( n, n, D, ldD, Dmat + nbl*nxn, n ) ; }

    void
    loadE( integer nbl, valueConstPointer E, integer ldE )
    { gecopy( n, n, E, ldE, Emat + nbl*nxn, n ) ; }

    t_Value & D( integer nbl, integer i, integer j )
    { return Dmat[ nbl*nxn + i + j*n] ; }

    t_Value const & D( integer nbl, integer i, integer j ) const
    { return Dmat[ nbl*nxn + i + j*n] ; }

    t_Value & E( integer nbl, integer i, integer j )
    { return Emat[ nbl*nxn + i + j*n] ; }

    t_Value const & E( integer nbl, integer i, integer j ) const
    { return Emat[ nbl*nxn + i + j*n] ; }

    // Border Bottom blocks
    void
    setZeroC()
    { zero( (nblock+1)*nxnb, Cmat, 1 ) ; }

    void
    loadC( integer nbl, valueConstPointer C, integer ldC )
    { gecopy( nb, n, C, ldC, Cmat + nbl*nxnb, nb ) ; }

    t_Value & C( integer nbl, integer i, integer j )
    { return Cmat[ nbl*nxnb + i + j*nb] ; }

    t_Value const & C( integer nbl, integer i, integer j ) const
    { return Cmat[ nbl*nxnb + i + j*nb] ; }

    // Border Right blocks
    void
    setZeroB()
    { zero( nblock*nxnb, Bmat, 1 ) ; }

    void
    loadB( integer nbl, valueConstPointer B, integer ldB )
    { gecopy( n, nb, B, ldB, Bmat + nbl*nxnb, n ) ; }

    t_Value & B( integer nbl, integer i, integer j )
    { return Bmat[ nbl*nxnb + i + j*n] ; }

    t_Value const & B( integer nbl, integer i, integer j ) const
    { return Bmat[ nbl*nxnb + i + j*n] ; }

    // Border RBblock
    void
    setZeroF()
    { zero( nb*nb, Fmat, 1 ) ; }

    void
    loadF( valueConstPointer F, integer ldF )
    { gecopy( nb, nb, F, ldF, Fmat, nb ) ; }

    t_Value & F( integer i, integer j )
    { return Fmat[ i + j*nb] ; }

    t_Value const & F( integer i, integer j ) const
    { return Fmat[ i + j*nb] ; }

    void
    loadBottom( valueConstPointer H0, integer ld0,
                valueConstPointer HN, integer ldN,
                valueConstPointer Hp, integer ldP,
                valueConstPointer Hq, integer ldQ ) ;

    t_Value & H( integer i, integer j )
    { return H0Npq[ i + j*(2*n+q)] ; }

    t_Value const & H( integer i, integer j ) const
    { return H0Npq[ i + j*(2*n+q)] ; }

    void
    factorize() ;

    void
    solve( valuePointer ) const ;

    void
    solve( integer      /* nrhs  */,
           valuePointer /* rhs   */,
           integer      /* ldRhs */ ) const ;

    // aux function
    void
    Mv( valueConstPointer x, valuePointer res ) const ;

    void
    dump_ccoord( ostream & stream ) const ;

  } ;
}

#endif
