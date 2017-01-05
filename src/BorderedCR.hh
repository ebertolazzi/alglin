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

#include "Alglin_threads.hh"
#ifndef BORDERED_CYCLIC_REDUCTION_MAX_THREAD
  #define BORDERED_CYCLIC_REDUCTION_MAX_THREAD 256
#endif

namespace alglin {

  //! available LU factorization code
  typedef enum {
    BORDERED_LAST_LU  = 0,
    BORDERED_LAST_QR  = 1,
    BORDERED_LAST_QRP = 2
  } BORDERED_LAST_Choice;

  typedef enum {
    BORDERED_LU  = 0,
    BORDERED_QR  = 1,
    BORDERED_QRP = 2
  } BORDERED_Choice;

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
  \*/

  /*\
   |
   |  Matrix structure
   |
   |                 n * (nblock+1)
   |    ___________________^____________________
   |   /                                        \
   |    n   n   n                              n  qx  nx
   |  +---+---+---+----.................-----+---+---+---+   -+
   |  | D | E |   |                          |   |   | B | n  |
   |  +---+---+---+                     -----+---+---+---+    |
   |  |   | D | E |                          |   |   | B | n  |
   |  +---+---+---+---+                 -----+---+---+---+    |
   |  |   |   | D | E |                      |   |   | B | n  |
   |  +---+---+---+---+                 -----+---+---+---+    |
   |  :                                                  :    |
   |  :                                                  :    |
   |  :                                                  :     > n * nblock
   |  :                                                  :    |
   |  :                                                  :    |
   |  :                              +---+---+---+---+---+    |
   |  :                              | D | E |   |   | B | n  |
   |  :                              +---+---+---+---+---+    |
   |  :                                  | D | E |   | B | n  |
   |  +---+---+---................---+---+---+---+---+---+   -+
   |  |   |   |                          |   |   |   |   |    |
   |  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
   |  |   |   |                          |   |   |   |   |    |
   |  +---+---+---................---+---+---+---+---+---+   -+
   |  | C | C |                      | C | C | C | Cq| F |    | nr
   |  +---+---+---................---+---+---+---+---+---+   -+
   |                                             nr*qx
  \*/

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
    integer qr, qx ; //!< extra BC
    integer nr, nx ; //!< border size

    // some derived constants
    integer n_x_2 ;
    integer n_x_n ;
    integer n_x_nx ;
    integer nr_x_n ;
    integer nr_x_nx ;
    integer Nr, Nc ;
    integer Tsize ;
    
    BORDERED_LAST_Choice last_selected ;
    BORDERED_Choice      selected ;

    void
    buildT( integer           nth,
            valueConstPointer TOP,
            valueConstPointer BOTTOM,
            valuePointer      T,
            integer *         iperm ) const ;

    void
    applyT( integer           nth,
            valueConstPointer T,
            integer const *   iperm,
            valuePointer      TOP,
            integer           ldTOP,
            valuePointer      BOTTOM,
            integer           ldBOTTOM,
            integer           ncol ) const ;

    void
    applyT( integer           nth,
            valueConstPointer T,
            integer const *   iperm,
            valuePointer      TOP,
            valuePointer      BOTTOM ) const ;

    /*
    //    __         _           _
    //   / _|__ _ __| |_ ___ _ _(_)______
    //  |  _/ _` / _|  _/ _ \ '_| |_ / -_)
    //  |_| \__,_\__|\__\___/_| |_/__\___|
    */

    void
    factorize_block( integer nth ) ;

    void
    factorize_reduced() ;
    
    /*
    //    __                            _
    //   / _|___ _ ___ __ ____ _ _ _ __| |
    //  |  _/ _ \ '_\ V  V / _` | '_/ _` |
    //  |_| \___/_|  \_/\_/\__,_|_| \__,_|
    */

    void
    forward( integer nth, valuePointer x ) const ;

    void
    forward_n( integer      nth,
               integer      nrhs,
               valuePointer rhs,
               integer      ldRhs ) const ;

    void
    forward_reduced( valuePointer x ) const ;

    void
    forward_n_reduced( integer      nrhs,
                       valuePointer rhs,
                       integer      ldRhs ) const ;

    /*
    //   _             _                       _
    //  | |__  __ _ __| |____ __ ____ _ _ _ __| |
    //  | '_ \/ _` / _| / /\ V  V / _` | '_/ _` |
    //  |_.__/\__,_\__|_\_\ \_/\_/\__,_|_| \__,_|
    */

    void
    backward( integer nth, valuePointer x ) const ;

    void
    backward_reduced( valuePointer x ) const ;

    void
    backward_n( integer      nth,
                integer      nrhs,
                valuePointer rhs,
                integer      ldRhs ) const ;

    void
    backward_n_reduced( integer      nrhs,
                        valuePointer rhs,
                        integer      ldRhs ) const ;

    void
    load_and_factorize_last() ;

    /*
    //   _         _
    //  | |__ _ __| |_
    //  | / _` (_-<  _|
    //  |_\__,_/__/\__|
    */

    void
    solve_last( valuePointer ) const ;

    void
    solve_last( integer      nrhs,
                valuePointer rhs,
                integer      ldRhs ) const ;

    valuePointer H0Nqp ;
    valuePointer Bmat, Cmat, Cqmat, Dmat, Emat, Fmat, WorkT, WorkQR ;
    
    // working block
    valuePointer Tmat, Ttau, Work ;
    integer      *Rperm, *Cperm, Lwork, LworkT, LworkQR ;

    // last block
    valuePointer Hmat, Htau ;
    integer      *Hperm, *Hswaps, *Tperm ;

    integer      *iBlock, *kBlock ;

    integer maxThread, usedThread ;
    mutable valuePointer xb_thread ;
    mutable std::thread threads[BORDERED_CYCLIC_REDUCTION_MAX_THREAD] ;
    mutable SpinLock spin ;

  public:

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    explicit
    BorderedCR( integer nth = integer(std::thread::hardware_concurrency()) )
    #else
    explicit
    BorderedCR()
    #endif
    : baseValue("BorderedCR_values")
    , baseInteger("BorderedCR_integers")
    , nblock(0)
    , n(0)
    , qr(0)
    , qx(0)
    , nr(0)
    , nx(0)
    , n_x_2(0)
    , n_x_n(0)
    , n_x_nx(0)
    , nr_x_n(0)
    , nr_x_nx(0)
    , Nr(0)
    , Nc(0)
    , last_selected(BORDERED_LAST_LU)
    , selected(BORDERED_LU)
    , H0Nqp(nullptr)
    , Bmat(nullptr)
    , Cmat(nullptr)
    , Cqmat(nullptr)
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
    {
      #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
      ALGLIN_ASSERT( nth > 0 && nth <= BORDERED_CYCLIC_REDUCTION_MAX_THREAD,
                     "Bad number of thread specification [" << nth << "]\n"
                     "must be a number > 0 and <= " << BORDERED_CYCLIC_REDUCTION_MAX_THREAD ) ;
      maxThread = nth ;
      #else
      maxThread = 1 ;
      #endif
    }

    virtual ~BorderedCR()
    {}

    //! load matrix in the class
    void
    allocate( integer _nblock,
              integer _n,
              integer _qr,
              integer _qx,
              integer _nr,
              integer _nx ) ;

    void select_LU()  { selected = BORDERED_LU ; }
    void select_QR()  { selected = BORDERED_QR ; }
    void select_QRP() { selected = BORDERED_QRP ; }

    void select_last_LU()  { last_selected = BORDERED_LAST_LU ; }
    void select_last_QR()  { last_selected = BORDERED_LAST_QR ; }
    void select_last_QRP() { last_selected = BORDERED_LAST_QRP ; }

    // filling bidiagonal part of the matrix
    // -------------------------------------------------------------------------
    void fillZero() ;

    void
    loadD( integer nbl, valueConstPointer D, integer ldD )
    { gecopy( n, n, D, ldD, Dmat + nbl*n_x_n, n ) ; }

    t_Value & D( integer nbl, integer i, integer j )
    { return Dmat[ nbl*n_x_n + i + j*n] ; }

    t_Value const & D( integer nbl, integer i, integer j ) const
    { return Dmat[ nbl*n_x_n + i + j*n] ; }

    void
    loadE( integer nbl, valueConstPointer E, integer ldE )
    { gecopy( n, n, E, ldE, Emat + nbl*n_x_n, n ) ; }

    t_Value & E( integer nbl, integer i, integer j )
    { return Emat[ nbl*n_x_n + i + j*n] ; }

    t_Value const & E( integer nbl, integer i, integer j ) const
    { return Emat[ nbl*n_x_n + i + j*n] ; }

    void
    loadDE( integer nbl, valueConstPointer DE, integer ldDE ) {
      gecopy( n, n, DE, ldDE, Dmat + nbl*n_x_n, n ) ; DE += n*ldDE ;
      gecopy( n, n, DE, ldDE, Emat + nbl*n_x_n, n ) ;
    }

    // -------------------------------------------------------------------------
    // Border Bottom blocks
    void
    setZeroC()
    { zero( (nblock+1)*nr_x_n, Cmat, 1 ) ; }

    void
    loadC( integer nbl, valueConstPointer C, integer ldC )
    { gecopy( nr, n, C, ldC, Cmat + nbl*nr_x_n, nr ) ; }

    void
    addtoC( integer nbl, valueConstPointer C, integer ldC ) {
      valuePointer CC = Cmat + nbl*nr_x_n ;
      geadd( nr, n, 1.0, C, ldC, 1.0, CC, nr, CC, nr ) ;
    }

    // add to block nbl and nbl+1
    void
    addtoC2( integer nbl, valueConstPointer C, integer ldC ) {
      valuePointer CC = Cmat + nbl*nr_x_n ;
      geadd( nr, n_x_2, 1.0, C, ldC, 1.0, CC, nr, CC, nr ) ;
    }

    t_Value & C( integer nbl, integer i, integer j )
    { return Cmat[ nbl*nr_x_n + i + j*nr ] ; }

    t_Value const & C( integer nbl, integer i, integer j ) const
    { return Cmat[ nbl*nr_x_n + i + j*nr ] ; }

    // -------------------------------------------------------------------------
    // Border Right blocks
    void
    setZeroB()
    { zero( nblock*n_x_nx, Bmat, 1 ) ; }

    void
    loadB( integer nbl, valueConstPointer B, integer ldB )
    { gecopy( n, nx, B, ldB, Bmat + nbl*n_x_nx, n ) ; }

    void
    addtoB( integer nbl, valueConstPointer B, integer ldB ) {
      valuePointer BB = Bmat + nbl*n_x_nx ;
      geadd( n, nx, 1.0, B, ldB, 1.0, BB, n, BB, n ) ;
    }

    t_Value & B( integer nbl, integer i, integer j )
    { return Bmat[ nbl*n_x_nx + i + j*n ] ; }

    t_Value const & B( integer nbl, integer i, integer j ) const
    { return Bmat[ nbl*n_x_nx + i + j*n ] ; }

    void
    loadDEB( integer nbl, valueConstPointer DEB, integer ldDEB ) {
      gecopy( n, n,  DEB, ldDEB, Dmat + nbl*n_x_n,  n  ) ; DEB += n*ldDEB ;
      gecopy( n, n,  DEB, ldDEB, Emat + nbl*n_x_n,  n  ) ; DEB += n*ldDEB ;
      gecopy( n, nx, DEB, ldDEB, Bmat + nbl*n_x_nx, nx ) ;
    }

    // -------------------------------------------------------------------------

    void
    setZeroF()
    { zero( nr*nx, Fmat, 1 ) ; }

    void
    loadF( valueConstPointer F, integer ldF )
    { gecopy( nr, nx, F, ldF, Fmat, nx ) ; }

    t_Value & F( integer i, integer j )
    { return Fmat[ i + j*nr ] ; }

    t_Value const & F( integer i, integer j ) const
    { return Fmat[ i + j*nr ] ; }

    // -------------------------------------------------------------------------

    void
    setZeroCq()
    { zero( nr*qx, Cqmat, 1 ) ; }

    void
    loadCq( valueConstPointer Cq, integer ldC )
    { gecopy( nr, qx, Cq, ldC, Cqmat, nr ) ; }

    t_Value & Cq( integer i, integer j )
    { return Cqmat[ i + j*nr ] ; }

    t_Value const & Cq( integer i, integer j ) const
    { return Cqmat[ i + j*nr ] ; }

    void
    loadCqF( valueConstPointer CqF, integer ldCF ) {
      gecopy( nr, qx, CqF, ldCF, Cqmat, nr ) ; CqF += qx*ldCF ;
      gecopy( nr, nx, CqF, ldCF, Fmat,  nr ) ;
    }

    // -------------------------------------------------------------------------

    void
    loadBottom( valueConstPointer H0, integer ld0,
                valueConstPointer HN, integer ldN,
                valueConstPointer Hq, integer ldQ,
                valueConstPointer Hp, integer ldP ) ;

    void
    loadBottom( valueConstPointer _H0Nqp, integer ldH ) {
      integer nq = n+qr ;
      gecopy( nq, Nc, _H0Nqp, ldH, H0Nqp, nq ) ;
    }

    t_Value & H( integer i, integer j )
    { return H0Nqp[ i + j*(n+qr) ] ; }

    t_Value const & H( integer i, integer j ) const
    { return H0Nqp[ i + j*(n+qr) ] ; }

    void
    factorize() ;

    void
    solve( valuePointer x ) const ;

    void
    solve( integer nrhs, valuePointer rhs, integer ldRhs ) const ;

    // aux function
    void
    Mv( valueConstPointer x, valuePointer res ) const ;

    void
    dump_ccoord( ostream & stream ) const ;

  } ;
}

#endif
