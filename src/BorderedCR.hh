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

  //! available LU factorization code
  typedef enum {
    BORDERED_LAST_LU  = 0,
    BORDERED_LAST_QR  = 1,
    BORDERED_LAST_QRP = 2
  } BORDERED_LAST_Choice;

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
    integer N ;
    
    BORDERED_LAST_Choice last_selected ;

    void
    buildT( valueConstPointer TOP,
            valueConstPointer BOTTOM,
            valuePointer      T,
            integer *         iperm ) const ;

    void
    applyT( valueConstPointer T,
            integer const *   iperm,
            valuePointer      TOP,
            integer           ldTOP,
            valuePointer      BOTTOM,
            integer           ldBOTTOM,
            integer           ncol ) const ;

    void
    applyT( valueConstPointer T,
            integer const *   iperm,
            valuePointer      TOP,
            valuePointer      BOTTOM ) const ;
    
    void
    load_last_block() ;

    void
    factorize_LU() ;

    void
    factorize_last_LU() ;

    void
    factorize_last_QR() ;

    void
    factorize_last_QRP() ;

    void
    solve_LU( valuePointer ) const ;

    void
    solve_LU( integer      /* nrhs  */,
              valuePointer /* rhs   */,
              integer      /* ldRhs */ ) const ;

    void
    solve_last_LU( valuePointer ) const ;

    void
    solve_last_QR( valuePointer ) const ;

    void
    solve_last_QRP( valuePointer ) const ;

    void
    solve_last_LU( integer      /* nrhs  */,
                   valuePointer /* rhs   */,
                   integer      /* ldRhs */ ) const ;

    void
    solve_last_QR( integer      /* nrhs  */,
                   valuePointer /* rhs   */,
                   integer      /* ldRhs */ ) const ;

    void
    solve_last_QRP( integer      /* nrhs  */,
                    valuePointer /* rhs   */,
                    integer      /* ldRhs */ ) const ;

    /*
    //
    //  Matrix structure
    //
    //                 n * (nblock+1)
    //    ___________________^____________________
    //   /                                        \
    //    n   n   n                              n   q  nb
    //  +---+---+---+----.................-----+---+---+---+   -+
    //  | D | E |   |                          |   |   | B | n  |
    //  +---+---+---+                     -----+---+---+---+    |
    //  |   | D | E |                          |   |   | B | n  |
    //  +---+---+---+---+                 -----+---+---+---+    |
    //  |   |   | D | E |                      |   |   | B | n  |
    //  +---+---+---+---+                 -----+---+---+---+    |
    //  :                                                  :    |
    //  :                                                  :    |
    //  :                                                  :     > n * nblock
    //  :                                                  :    |
    //  :                                                  :    |
    //  :                              +---+---+---+---+---+    |
    //  :                              | D | E |   |   | B | n  |
    //  :                              +---+---+---+---+---+    |
    //  :                                  | D | E |   | B | n  |
    //  +---+---+---................---+---+---+---+---+---+   -+
    //  |   |   |                          |   |   |   |   |    |
    //  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+q
    //  |   |   |                          |   |   |   |   |    |
    //  +---+---+---................---+---+---+---+---+---+   -+
    //  | C | C |                      | C | C | C | 0 | F |    | nb
    //  +---+---+---................---+---+---+---+---+---+   -+
    */

    valuePointer H0Nqp ;
    valuePointer Bmat, Cmat, Dmat, Emat, Fmat ;
    
    // working block
    valuePointer Tmat, Ttau, Work ;
    integer      *Rperm, *Cperm, Lwork ;

    // last block
    valuePointer Hmat, Htau ;
    integer      *Hperm, *Hswaps ;

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
    , N(0)
    , last_selected(BORDERED_LAST_QRP)
    , H0Nqp(nullptr)
    , Bmat(nullptr)
    , Cmat(nullptr)
    , Dmat(nullptr)
    , Emat(nullptr)
    , Fmat(nullptr)
    , Tmat(nullptr)
    , Ttau(nullptr)
    , Work(nullptr)
    , Rperm(nullptr)
    , Cperm(nullptr)
    , Lwork(0)
    , Hmat(nullptr)
    , Htau(nullptr)
    , Hperm(nullptr)
    {}

    virtual ~BorderedCR()
    {}

    //! load matrix in the class
    void
    allocate( integer _nblock,
              integer _n,
              integer _q,
              integer _nb ) ;

    void select_last_LU()  { last_selected = BORDERED_LAST_LU ; }
    void select_last_QR()  { last_selected = BORDERED_LAST_QR ; }
    void select_last_QRP() { last_selected = BORDERED_LAST_QRP ; }

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
                valueConstPointer Hq, integer ldQ,
                valueConstPointer Hp, integer ldP ) ;

    void
    loadBottom( valueConstPointer _H0Nqp, integer ldH ) {
      integer nq = n+q ;
      gecopy( nq, N, _H0Nqp, ldH, H0Nqp, nq ) ;
    }

    t_Value & H( integer i, integer j )
    { return H0Nqp[ i + j*(n+q)] ; }

    t_Value const & H( integer i, integer j ) const
    { return H0Nqp[ i + j*(n+q)] ; }

    void
    factorize() {
      factorize_LU() ;
      load_last_block() ;
      switch ( last_selected ) {
        case BORDERED_LAST_LU:  factorize_last_LU()  ; break ;
        case BORDERED_LAST_QR:  factorize_last_QR()  ; break ;
        case BORDERED_LAST_QRP: factorize_last_QRP() ; break ;
      }
    }

    void
    solve( valuePointer x ) const
    { solve_LU( x ) ; }

    void
    solve( integer      nrhs,
           valuePointer rhs,
           integer      ldRhs ) const
    { solve_LU( nrhs, rhs, ldRhs ) ; }

    // aux function
    void
    Mv( valueConstPointer x, valuePointer res ) const ;

    void
    dump_ccoord( ostream & stream ) const ;

  } ;
}

#endif
