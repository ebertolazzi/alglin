/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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

#include "Alglin.hh"
#include "Alglin_Eigen.hh"

#define BABD_LAST_ERROR( ... ) m_last_error = fmt::format(__VA_ARGS__ )

#define BABD_LAST_ERROR_LOCK( ... ) m_spin.lock(); m_last_error = fmt::format(__VA_ARGS__ ); m_spin.unlock()

namespace alglin {

  using std::exception;
  using std::min;
  using std::max;
  using std::swap;

  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::allocate(
    integer nblock,
    integer n,
    integer qr,
    integer qx,
    integer nr,
    integer nx
  ) {

    if ( m_number_of_blocks == nblock &&
         m_block_size       == n      &&
         m_qr               == qr     &&
         m_qx               == qx     &&
         m_nr               == nr     &&
         m_nx               == nx ) return;

    try {

      m_number_of_blocks = nblock;
      m_block_size       = n;

      m_qr    = qr;
      m_qx    = qx;
      m_nr    = nr;
      m_nx    = nx;
      n_x_2   = n*2;
      n_x_n   = n*n;
      nr_x_n  = n*m_nr;
      n_x_nx  = n*m_nx;
      nr_x_nx = m_nr*m_nx;
      nr_x_qx = m_nr*m_qx;
      m_Nr    = n_x_2+m_nr+m_qr;
      m_Nc    = n_x_2+m_nx+m_qx;
      m_Tsize = 2*n_x_n+n;

      integer NBLK = max(1,m_max_parallel_block);

      integer N = max(m_Nr,m_Nc);
      m_Work_Lapack_size = max(N,2*max(n_x_n,max(nr_x_n,n_x_nx)));

      real_type tmp; // get optimal allocation
      integer info = alglin::geqrf( m_Nr, m_Nc, nullptr, m_Nr, nullptr, &tmp, -1 );
      UTILS_ASSERT(
        info == 0,
        "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
      );
      if ( m_Work_Lapack_size < integer(tmp) ) m_Work_Lapack_size = integer(tmp);

      info = geqp3( m_Nr, m_Nc, nullptr, m_Nr, nullptr, nullptr, &tmp, -1 );
      UTILS_ASSERT(
        info == 0,
        "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
      );
      if ( m_Work_Lapack_size < integer(tmp) ) m_Work_Lapack_size = integer(tmp);

      info = alglin::geqrf( n_x_2, n, nullptr, n_x_2, nullptr, &tmp, -1 );
      UTILS_ASSERT(
        info == 0,
        "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
      );
      if ( m_Work_Lapack_size < integer(tmp) ) m_Work_Lapack_size = integer(tmp);

      info = alglin::geqp3( n_x_2, n, nullptr, n_x_2, nullptr, nullptr, &tmp, -1 );
      UTILS_ASSERT(
        info == 0,
        "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
      );
      if ( m_Work_Lapack_size < integer(tmp) ) m_Work_Lapack_size = integer(tmp);

      integer nnz = (n*(nr+2*n+nx+1)+m_Tsize)*nblock +
                    (n + m_Nr + qr)*m_Nc + (n + qx)*nr +
                    (m_Work_Lapack_size+nr+nr_x_nx)*NBLK;
      integer innz = nblock*n + (3+n)*NBLK;

      m_mem.reallocate( size_t(nnz) );
      m_mem_int.reallocate( size_t(innz) );

      m_Bmat  = m_mem( size_t(nblock*n_x_nx) );
      m_Cmat  = m_mem( size_t((nblock+1)*nr_x_n) );
      m_Cqmat = m_mem( size_t(nr_x_qx) );
      m_Dmat  = m_mem( size_t(nblock*n_x_n) );
      m_Emat  = m_mem( size_t(nblock*n_x_n) );
      m_Fmat  = m_mem( size_t(NBLK*nr_x_nx) );

      m_Tmat  = m_mem( size_t(nblock*m_Tsize) );
      m_Ttau  = m_mem( size_t(nblock*n) );
      m_Hmat  = m_mem( size_t(m_Nr*m_Nc) );
      m_H0Nqp = m_mem( size_t((n+qr)*m_Nc) );

      m_Perm   = m_mem_int( size_t(nblock*n) );
      m_iBlock = m_mem_int( size_t(2*NBLK) );
      m_kBlock = m_mem_int( size_t(NBLK) );

      m_perm_thread.resize( size_t(NBLK) );
      m_xb_thread.resize( size_t(NBLK) );
      m_Work_Lapack_thread.resize( size_t(NBLK) );
      m_Work_T_thread.resize( size_t(NBLK) );

      // precompute partition for parallel computation
      for ( size_t nt = 0; nt < size_t(NBLK); ++nt ) {
        m_perm_thread[nt]        = m_mem_int( size_t(n) );
        m_xb_thread[nt]          = m_mem( size_t(nr) );
        m_Work_Lapack_thread[nt] = m_mem( size_t(m_Work_Lapack_size) );
        m_Work_T_thread[nt].resize(100);
      }

      // calcolo partizionamento blocchi in modo che ogni thread
      // abbia almeno 2 blocchi di righe da processare.
      m_used_parallel_block = NBLK;
      m_reduced_nblk        = 2*NBLK;
      while ( m_reduced_nblk > m_number_of_blocks )
        { m_reduced_nblk -= 2; --m_used_parallel_block; }

      if ( m_used_parallel_block <= 1 ) {
        m_used_parallel_block = 1;
        m_iBlock[0]           = 0;
        m_iBlock[1]           = m_number_of_blocks;
      } else {
        m_iBlock[0] = 0;
        m_iBlock[1] = static_cast<integer>(m_number_of_blocks/m_used_parallel_block);
        for ( integer nt = 1; nt < m_used_parallel_block; ++nt ) {
          m_iBlock[2*nt+0] = m_iBlock[2*nt-1]+1;
          m_iBlock[2*nt+1] = static_cast<integer>(((nt+1)*m_number_of_blocks)/m_used_parallel_block);
        }
        --m_reduced_nblk;
      }

      m_mem.must_be_empty("BorderedCR::allocate, values");
      m_mem_int.must_be_empty("BorderedCR::allocate, integer");

    } catch ( exception const & err ) {
      m_last_error = err.what();
      throw;
    } catch (...) {
      m_last_error = "BorderedCR::allocate, unknown error\n";
      throw;
    }
  }

  /*\
   |       _
   |    __| |_   _ _ __
   |   / _` | | | | '_ \
   |  | (_| | |_| | |_) |
   |   \__,_|\__,_| .__/
   |              |_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::dup( BorderedCR const & M ) {
    m_max_parallel_block = M.m_max_parallel_block;
    allocate( M.m_number_of_blocks, M.m_block_size, M.m_qr, M.m_qx, M.m_nr, M.m_nx );

    integer const & nblk {m_number_of_blocks};
    integer const & n    {m_block_size};

    Copy_n( M.m_Bmat,  nblk*n_x_nx,     m_Bmat  );
    Copy_n( M.m_Cmat,  (nblk+1)*nr_x_n, m_Cmat  );
    Copy_n( M.m_Dmat,  nblk*n_x_n,      m_Dmat  );
    Copy_n( M.m_Emat,  nblk*n_x_n,      m_Emat  );
    Copy_n( M.m_Fmat,  nr_x_nx,         m_Fmat  );
    Copy_n( M.m_H0Nqp, (n+m_qr)*m_Nc,   m_H0Nqp );
    Copy_n( M.m_Cqmat, nr_x_qx,         m_Cqmat );
    m_TP = &m_TP_fake;
  }

  /*\
   |   _        __
   |  (_)_ __  / _| ___
   |  | | '_ \| |_ / _ \
   |  | | | | |  _| (_) |
   |  |_|_| |_|_|  \___/
  \*/

  template <typename t_Value>
  string
  BorderedCR<t_Value>::info( string_view indent ) const {
    return fmt::format(
      "{0}rows   = {1}\n"
      "{0}cols   = {2}\n"
      "{0}nblock = {3}\n"
      "{0}n      = {4}\n"
      "{0}qr     = {5}\n"
      "{0}qx     = {6}\n"
      "{0}nr     = {7}\n"
      "{0}nx     = {8}\n"
      "{0}treads = {9}/nblk:{10}/factorize:{11}/solve:{12}\n"
      "{0}intern = {13}\n"
      "{0}last   = {14}\n"
      "{0}last2  = {14}\n",
      indent,
      nrows(), ncols(), m_number_of_blocks,
      m_block_size, m_qr, m_qx, m_nr, m_nx,
      m_TP->thread_count(),
      m_used_parallel_block,
      m_factorize_use_thread,
      m_solve_use_thread,
      choice_to_string( m_selected ),
      choice_to_string( m_last_block_selected ),
      choice_to_string( m_last_block_selected2 )
    );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_D() {
    Zero_n( m_Dmat, m_number_of_blocks*n_x_n );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_E() {
    Zero_n( m_Emat, m_number_of_blocks*n_x_n );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_B() {
    Zero_n( m_Bmat, m_number_of_blocks*n_x_nx );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_F() {
    Zero_n( m_Fmat, nr_x_nx );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_H() {
    integer const & n{m_block_size};
    Zero_n( m_H0Nqp, (n+m_qr)*m_Nc );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_C() {
    integer const & nblk{m_number_of_blocks};
    Zero_n( m_Cmat, (nblk+1)*nr_x_n );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zero_Cq() {
    Zero_n( m_Cqmat, nr_x_qx );
  }

  /*\
   |   _                 _ ____        _   _
   |  | | ___   __ _  __| | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` |  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_|____/ \___/ \__|\__\___/|_| |_| |_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom(
    real_type const H0[], integer ld0,
    real_type const HN[], integer ldN,
    real_type const Hq[], integer ldQ,
    real_type const Hp[], integer ldP
  ) {
    integer const & n{m_block_size};
    // (n+qr) x ( n + n + qx + nx )
    integer     m{n + m_qr};
    real_type * H{m_H0Nqp};
    GEcopy( m, n,    H0, ld0, H, m ); H += m * n;
    GEcopy( m, n,    HN, ldN, H, m ); H += m * n;
    GEcopy( m, m_qx, Hq, ldQ, H, m ); H += m * m_qx;
    GEcopy( m, m_nx, Hp, ldP, H, m );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom(
    MatW const & H0,
    MatW const & HN,
    MatW const & Hq,
    MatW const & Hp
  ) {
    integer const & n{m_block_size};
    integer m{n + m_qr};

    bool ok{true};
    if ( ! (H0.nrows() == m && H0.ncols() == n) ) {
      ok = false;
      BABD_LAST_ERROR(
        "BorderedCR::load_bottom, bad dimension size(H0) = {} x {} expected {} x {}\n",
        H0.nrows(), H0.ncols(), m, n
      );
    } else if ( ! (HN.nrows() == m && HN.ncols() == n) ) {
      ok = false;
      BABD_LAST_ERROR(
        "BorderedCR::load_bottom, bad dimension size(HN) = {} x {} expected {} x {}\n",
        HN.nrows(), HN.ncols(), m, n
      );
    } else if ( ! (Hq.nrows() == m && Hq.ncols() == m_qx) ) {
      ok = false;
      BABD_LAST_ERROR(
        "BorderedCR::load_bottom, bad dimension size(Hq) = {} x {} expected {} x {}\n",
        Hq.nrows(), Hq.ncols(), m, m_qx
      );
    } else if ( ! (Hp.nrows() == m && Hp.ncols() == m_nx) ) {
      ok = false;
      BABD_LAST_ERROR(
        "BorderedCR::load_bottom, bad dimension size(Hp) = {} x {} expected {} x {}\n",
        Hp.nrows(), Hp.ncols(), m, m_nx
      );
    } else {
      load_bottom(
        H0.data(), H0.ldim(),
        HN.data(), HN.ldim(),
        Hq.data(), Hq.ldim(),
        Hp.data(), Hp.ldim()
      );
    }
    if ( !ok ) Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom( real_type const H0Nqp[], integer ldH ) {
    integer const & n{m_block_size};
    integer nq{n + m_qr};
    GEcopy( nq, m_Nc, H0Nqp, ldH, m_H0Nqp, nq );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom( MatW const & H ) {
    integer const & n{m_block_size};
    integer m{n + m_qr};
    if ( H.nrows() == m && H.ncols() == m_Nc ) {
      load_bottom( H.data(), H.ldim() );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_bottom, bad dimension size(H) = {} x {} expected {} x {}\n",
        H.nrows(), H.ncols(), m, m_Nc
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom2(
    real_type const C0[], integer ld0,
    real_type const CN[], integer ldN,
    real_type const Cq[], integer ldCq,
    real_type const F[],  integer ldF
  ) {
    integer const & nblk {m_number_of_blocks};
    integer const & n    {m_block_size};
    // (n+qr) x ( n + n + qx + nx )
    GEcopy( m_nr, n,    C0, ld0,  m_Cmat,            m_nr );
    GEcopy( m_nr, n,    CN, ldN,  m_Cmat+nblk*n_x_n, m_nr );
    GEcopy( m_nr, m_qx, Cq, ldCq, m_Cqmat,           m_nr );
    GEcopy( m_nr, m_nx, F,  ldF,  m_Fmat,            m_nr );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom2(
    MatW const & C0,
    MatW const & CN,
    MatW const & Cq,
    MatW const & F
  ) {
    integer const & n{m_block_size};

    bool ok{ C0.nrows() == m_nr && C0.ncols() == n };
    if ( !ok ) {
      BABD_LAST_ERROR(
        "BorderedCR::load_bottom2, bad dimension size(C0) = {} x {} expected {} x {}\n",
        C0.nrows(), C0.ncols(), m_nr, n
      );
    } else {
      ok = CN.nrows() == m_nr && CN.ncols() == n;
      if ( !ok ) {
        BABD_LAST_ERROR(
          "BorderedCR::load_bottom2, bad dimension size(CN) = {} x {} expected {} x {}\n",
          CN.nrows(), CN.ncols(), m_nr, n
        );
      } else {
        ok = Cq.nrows() == m_nr && Cq.ncols() == m_qx;
        if ( !ok ) {
          BABD_LAST_ERROR(
            "BorderedCR::load_bottom2, bad dimension size(Cq) = {} x {} expected {} x {}\n",
            Cq.nrows(), Cq.ncols(), m_nr, m_qx
          );
        } else {
          ok = F.nrows() == m_nr && F.ncols() == m_nx;
          if ( !ok ) {
            BABD_LAST_ERROR(
              "BorderedCR::load_bottom2, bad dimension size(F) = {} x {} expected {} x {}\n",
              F.nrows(), F.ncols(), m_nr, m_nx
            );
          }
        }
      }
    }
    if ( ok ) {
      load_bottom2(
        C0.data(), C0.ldim(),
        CN.data(), CN.ldim(),
        Cq.data(), Cq.ldim(),
        F.data(),  F.ldim()
      );
    } else {
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_bottom2( MatW const & H ) {
    integer   const & nblock{m_number_of_blocks};
    real_type const * ptr{H.data()};
    integer           ld{H.ldim()};
    this->load_C(  0,      ptr, ld ); ptr += nr_x_n;
    this->load_C(  nblock, ptr, ld ); ptr += nr_x_n;
    this->load_Cq( ptr, ld ); ptr += m_nr * m_qx;
    this->load_F(  ptr, ld );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_B( integer nbl, real_type const B[], integer ldB ) {
    integer const & n{m_block_size};
    GEcopy( n, m_nx, B, ldB, m_Bmat + nbl*n_x_nx, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_B( integer nbl, MatW const & B ) {
    integer const & n{m_block_size};
    if ( B.nrows() == n && B.ncols() == m_nx ) {
      GEcopy( n, m_nx, B.data(), B.ldim(), m_Bmat + nbl*n_x_nx, n );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_B( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
        nbl, B.nrows(), B.ncols(), n, m_nx
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_B( integer nbl, real_type const B[], integer ldB ) {
    integer const & n{m_block_size};
    real_type * BB{m_Bmat + nbl*n_x_nx};
    GEadd( n, m_nx, B, ldB, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_B( integer nbl, MatW const & B ) {
    integer const & n{m_block_size};
    if ( B.nrows() == n && B.ncols() == m_nx ) {
      real_type * BB{m_Bmat + nbl*n_x_nx};
      GEadd( n, m_nx, B.data(), B.ldim(), BB, n );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::add_to_B( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
        nbl, B.nrows(), B.ncols(), n, m_nx
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_C( integer nbl, real_type const C[], integer ldC ) {
    integer const & n{m_block_size};
    GEcopy( m_nr, n, C, ldC, m_Cmat + nbl*nr_x_n, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_C( integer nbl, MatW const & C ) {
    integer const & n{m_block_size};
    if ( C.nrows() == m_nr && C.ncols() == n ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEcopy( m_nr, n, C.data(), C.ldim(), CC, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_C( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
        nbl, C.nrows(), C.ncols(), m_nr, n
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_C( integer nbl, real_type const C[], integer ldC ) {
    integer const & n{m_block_size};
    if ( ldC >= m_nr ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEadd( m_nr, n, C, ldC, CC, m_nr );
    } else {
      BABD_LAST_ERROR( "BorderedCR::add_to_C( {}, C, ldC = {} ) bad ldC\n", nbl, ldC );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_C( integer nbl, MatW const & C ) {
    integer const & n{m_block_size};
    if ( C.nrows() == m_nr && C.ncols() == n ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEadd( m_nr, n, C.data(), C.ldim(), CC, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::add_to_C( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
        nbl, C.nrows(), C.ncols(), m_nr, n
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_C2( integer nbl, real_type const C[], integer ldC ) {
    if ( ldC >= m_nr ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEadd( m_nr, n_x_2, C, ldC, CC, m_nr );
    } else {
      BABD_LAST_ERROR( "BorderedCR::add_to_C2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_C2( integer nbl, MatW const & C ) {
    if ( C.nrows() == m_nr && C.ncols() == n_x_2 ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEadd( m_nr, n_x_2, C.data(), C.ldim(), CC, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::add_to_C2( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
        nbl, C.nrows(), C.ncols(), m_nr, n_x_2
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_D( integer nbl, real_type const D[], integer ldD ) {
    integer const & n{m_block_size};
    GEcopy( n, n, D, ldD, m_Dmat + nbl*n_x_n, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_D( integer nbl, MatW const & D ) {
    integer const & n{m_block_size};
    if ( D.nrows() == n && D.ncols() == n ) {
      real_type * DD{m_Dmat + nbl*n_x_n};
      GEcopy( n, n, D.data(), D.ldim(), DD, n );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_D( {}, D) bad dimension size(D) = {} x {} expected {} x {}\n",
        nbl, D.nrows(), D.ncols(), n, n
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_E( integer nbl, real_type const E[], integer ldE ) {
    integer const & n{m_block_size};
    real_type * EE{m_Emat + nbl*n_x_n};
    GEcopy( n, n, E, ldE, EE, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_E( integer nbl, MatW const & E ) {
    integer const & n{m_block_size};
    if ( E.nrows() == n && E.ncols() == n ) {
      real_type * EE{m_Emat + nbl*n_x_n};
      GEcopy( n, n, E.data(), E.ldim(), EE, n );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_E( {}, E) bad dimension size(E) = {} x {} expected {} x {}\n",
        nbl, E.nrows(), E.ncols(), n, n
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_DE( integer nbl, real_type const DE[], integer ldDE ) {
    integer const & n{m_block_size};
    real_type * DD{m_Dmat + nbl*n_x_n};
    real_type * EE{m_Emat + nbl*n_x_n};
    GEcopy( n, n, DE, ldDE, DD, n ); DE += n*ldDE;
    GEcopy( n, n, DE, ldDE, EE, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_DEB( integer nbl, real_type const DEB[], integer ldDEB ) {
    integer const & n{m_block_size};
    real_type * DD{m_Dmat + nbl*n_x_n};
    real_type * EE{m_Emat + nbl*n_x_n};
    real_type * BB{m_Bmat + nbl*n_x_nx};
    GEcopy( n, n,    DEB, ldDEB, DD, n ); DEB += n*ldDEB;
    GEcopy( n, n,    DEB, ldDEB, EE, n ); DEB += n*ldDEB;
    GEcopy( n, m_nx, DEB, ldDEB, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_F( real_type const F[], integer ldF ) {
    GEcopy( m_nr, m_nx, F, ldF, m_Fmat, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_F( MatW const & F ) {
    if ( F.nrows() == m_nr && F.ncols() == m_nx ) {
      GEcopy( m_nr, m_nx, F.data(), F.ldim(), m_Fmat, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_F(F) bad dimension size(F) = {} x {} expected {} x {}\n",
        F.nrows(), F.ncols(), m_nr, m_nx
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_F( real_type const F[], integer ldF ) {
    GEadd( m_nr, m_nx, F, ldF, m_Fmat, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_F( MatW const & F ) {
    if ( F.nrows() == m_nr && F.ncols() == m_nx ) {
      GEadd( m_nr, m_nx, F.data(), F.ldim(), m_Fmat, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::add_to_F(F) bad dimension size(F) = {} x {} expected {} x {}\n",
        F.nrows(), F.ncols(), m_nr, m_nx
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_Cq( real_type const Cq[], integer ldC ) {
    GEcopy( m_nr, m_qx, Cq, ldC, m_Cqmat, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_Cq( MatW const & Cq ) {
    if ( Cq.nrows() == m_nr && Cq.ncols() == m_qx ) {
      GEcopy( m_nr, m_qx, Cq.data(), Cq.ldim(), m_Cqmat, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::load_Cq(Cq) bad dimension size(Cq) = {} x {} expected {} x {}\n",
        Cq.nrows(), Cq.ncols(), m_nr, m_qx
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_CqF( real_type const CqF[], integer ldCF ) {
    integer const & nr{m_nr};
    integer const & nx{m_nx};
    integer const & qx{m_qx};
    GEcopy( nr, qx, CqF, ldCF, m_Cqmat, nr ); CqF += qx*ldCF;
    GEcopy( nr, nx, CqF, ldCF, m_Fmat,  nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_C2F( integer nbl, real_type const C2F[], integer ldC ) {
    if ( ldC >= m_nr ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEadd( m_nr, n_x_2, C2F,               ldC, CC,     m_nr );
      GEadd( m_nr, m_nx,  C2F + ldC * n_x_2, ldC, m_Fmat, m_nr );
    } else {
      BABD_LAST_ERROR( "BorderedCR::add_to_C2F( {}, C2F, ldC = {} ) bad ldC\n", nbl, ldC );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_to_C2F( integer nbl, MatW const & C2F ) {
    if ( C2F.nrows() == m_nr && C2F.ncols() == n_x_2+m_nx ) {
      real_type * CC{m_Cmat + nbl*nr_x_n};
      GEadd( m_nr, n_x_2, C2F.data(),                m_nr, CC,     m_nr );
      GEadd( m_nr, m_nx,  C2F.data() + m_nr * n_x_2, m_nr, m_Fmat, m_nr );
    } else {
      BABD_LAST_ERROR(
        "BorderedCR::add_to_C2F( {}, C2F) bad dimension size(C2F) = {} x {} expected {} x {}\n",
        nbl, C2F.nrows(), C2F.ncols(), m_nr, n_x_2+m_nx
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  /*\
   |   _____          _             _
   |  |  ___|_ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::factorize() {

    if ( m_debug ) check_matrix();

    m_matrix_is_factorized = false;
    integer const & nblock {m_number_of_blocks};
    integer const & npb    {m_used_parallel_block};
    m_ok_thread = true;
    if ( npb > 1 ) {
      Zero_n( m_Fmat + nr_x_nx, (npb-1) * nr_x_nx );
      if ( m_factorize_use_thread ) {
        for ( integer n_thread = 0; n_thread < npb; ++n_thread )
          m_TP->run( &BorderedCR<t_Value>::factorize_block, this, n_thread );
        m_TP->wait();
      } else {
        for ( integer n_thread = 0; n_thread < npb; ++n_thread )
          factorize_block( n_thread );
      }
      if ( nr_x_nx > 0 ) {
        MapMatrix<t_Value> FMAT( m_Fmat, m_nr, m_nx );
        for ( integer n_thread = 1; n_thread < npb; ++n_thread ) {
          MapMatrix<t_Value> WORK( m_Fmat + n_thread*nr_x_nx, m_nr, m_nx );
          FMAT.noalias() += WORK;
        }
      }
      if ( m_ok_thread ) m_ok_thread = factorize_reduced();
    } else {
      m_iBlock[0] = 0;
      m_iBlock[1] = nblock;
      factorize_block( 0 );
    }
    if ( m_ok_thread ) m_ok_thread = load_and_factorize_last();
    m_matrix_is_factorized = m_ok_thread;
    return m_ok_thread;
  }

  /*\
    Compute   / L \ (U) = / LU \ = / TOP    \
              \ M /       \ MU /   \ BOTTOM /
  \*/
  template <typename t_Value>
  bool
  BorderedCR<t_Value>::buildT(
    integer   n_thread,
    real_type T[],
    integer   iperm[]
  ) const {
    integer const & n{m_block_size};
    integer info{0};
    switch ( m_selected ) {
    case BORDERED_Choice::LU:
      info = alglin::getrf( n_x_2, n, T, n_x_2, iperm );
      if ( info != 0 ) {
        BABD_LAST_ERROR(
          "BorderedCR::buildT [LU] getrf(M={},N={},A,LDA={},iperm)::INFO = {}\n",
          n_x_2, n, n_x_2, info
        );
        return false;
      }
      break;
    case BORDERED_Choice::QR:
      info = alglin::geqrf(
        n_x_2, n, T, n_x_2, T+2*n_x_n,
        m_Work_Lapack_thread[size_t(n_thread)], m_Work_Lapack_size
      );
      if ( info != 0 ) {
        BABD_LAST_ERROR(
          "BorderedCR::buildT [QR] geqrf(M={},N={},A,LDA={},...)::INFO = {}\n",
          n_x_2, n, n_x_2, info
        );
        return false;
      }
      break;
    case BORDERED_Choice::QRP:
      {
        integer * P{m_perm_thread[size_t(n_thread)]};
        Fill_n( P, n, 0 );
        info = alglin::geqp3(
          n_x_2, n, T, n_x_2, P, T+2*n_x_n,
          m_Work_Lapack_thread[size_t(n_thread)], m_Work_Lapack_size
        );
        if ( info != 0 ) {
          BABD_LAST_ERROR(
            "BorderedCR::buildT [QRP] geqp3(M={},N={},A,LDA={},...)::INFO = {}\n",
            n_x_2, n, n_x_2, info
          );
          return false;
        }
        permutation_to_exchange( n, P, iperm );
      }
      break;
    }
    return true;
  }

  /*\
   APPLY
      / -M  I \ / L^(-1)   \ / TOP    \
      \  I  0 / \        I / \ BOTTOM /
      = / -M  I \ / L^(-1) TOP \
        \  I  0 / \  BOTTOM    /
      = / BOTTOM - M L^(-1) TOP \
        \ L^(-1) TOP            /

   ON
      / -M  I \ / L^(-1)   \ / L \ (U)
      \  I  0 / \        I / \ M /
        = / -M  I \ / I \ (U)
          \  I  0 / \ M /
        = / 0 \ (U)   / 0 \
          \ I /     = \ U /
  \*/
  template <typename t_Value>
  bool
  BorderedCR<t_Value>::applyT(
    integer         n_thread,
    real_type const T[],
    integer   const iperm[],
    real_type       TOP[],
    integer         ldTOP,
    real_type       BOTTOM[],
    integer         ldBOTTOM,
    integer         ncol
  ) const {
    integer const & n{m_block_size};
    real_type     * W{Work_T(n_thread,n_x_2*ncol)};

    GEcopy( n, ncol, TOP,    ldTOP,    W,   n_x_2 );
    GEcopy( n, ncol, BOTTOM, ldBOTTOM, W+n, n_x_2 );

    // Apply row interchanges to the right hand sides.
    integer info{0};
    switch ( m_selected ) {
    case BORDERED_Choice::LU:
      info = alglin::swaps( ncol, W, n_x_2, 0, n-1, iperm, 1 );
      if ( info != 0 ) {
        BABD_LAST_ERROR_LOCK( "BorderedCR::applyT INFO = {}\n", info );
        return false;
      }
      alglin::trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, ncol, 1.0, T, n_x_2, W, n_x_2
      );
      alglin::gemm(
        Transposition::NO,
        Transposition::NO,
        n, ncol, n,
        -1.0, T+n, n_x_2,  // M
              W,   n_x_2,  // L^(-1) TOP
         1.0, W+n, n_x_2   // TOP = BOTTOM - M L^(-1) TOP
      );
      break;
    case BORDERED_Choice::QR:
    case BORDERED_Choice::QRP:
      info = alglin::ormqr(
        SideMultiply::LEFT,
        Transposition::YES,
        n_x_2, ncol, // righe x colonne
        n,           // numero riflettori usati nel prodotto Q
        T, n_x_2,
        T+2*n_x_n,
        W, n_x_2,
        m_Work_Lapack_thread[size_t(n_thread)], m_Work_Lapack_size
      );
      break;
    }

    GEcopy( n, ncol, W+n, n_x_2, TOP,    ldTOP    );
    GEcopy( n, ncol, W,   n_x_2, BOTTOM, ldBOTTOM );
    return true;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::applyT(
    integer         n_thread,
    real_type const T[],
    integer   const iperm[],
    real_type       TOP[],
    real_type       BOTTOM[]
  ) const {

    integer const & n{m_block_size};
    real_type     * W{Work_T(n_thread,n_x_2)};
    Copy_n( TOP,    n, W   );
    Copy_n( BOTTOM, n, W+n );
    integer info{0};
    switch ( m_selected ) {
    case BORDERED_Choice::LU:
      // Apply row interchanges to the right hand sides.
      info = alglin::swaps( 1, W, n_x_2, 0, n-1, iperm, 1 );
      if ( info != 0 ) {
        BABD_LAST_ERROR_LOCK( "BorderedCR::applyT, swaps return INFO = {}\n", info );
        return false;
      }
      alglin::trsv(
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, T, n_x_2, W, 1
      );
      alglin::gemv(
        Transposition::NO,
        n, n,
        -1.0, T+n, n_x_2,
              W,       1,
         1.0, W+n, 1
      );
      break;
    case BORDERED_Choice::QR:
    case BORDERED_Choice::QRP:
      info = alglin::ormqr(
        SideMultiply::LEFT,
        Transposition::YES,
        n_x_2, 1, // righe x colonne
        n,        // numero riflettori usati nel prodotto Q
        T, n_x_2,
        T+2*n_x_n,
        W, n_x_2,
        m_Work_Lapack_thread[size_t(n_thread)], m_Work_Lapack_size
      );
      if ( info != 0 ) {
        BABD_LAST_ERROR_LOCK( "BorderedCR::applyT, ormqr return INFO = {}\n", info );
        return false;
      }
      break;
    }
    Copy_n( W+n, n, TOP   );
    Copy_n( W,   n, BOTTOM );
    return true;
  }

  /*\
   |    __            _             _             _     _            _
   |   / _| __ _  ___| |_ ___  _ __(_)_______    | |__ | | ___   ___| | __
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \   | '_ \| |/ _ \ / __| |/ /
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/   | |_) | | (_) | (__|   <
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|___|_.__/|_|\___/ \___|_|\_\
   |                                        |_____|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_block( integer n_thread ) {

    integer const & n{m_block_size};

    integer iblock { m_iBlock[2*n_thread+0] };
    integer eblock { m_iBlock[2*n_thread+1] };
    integer nblk   { eblock - iblock };

    real_type * Bmat0 { m_Bmat + iblock   * n_x_nx  };
    real_type * Cmat0 { m_Cmat + iblock   * nr_x_n  };
    real_type * Dmat0 { m_Dmat + iblock   * n_x_n   };
    real_type * Emat0 { m_Emat + iblock   * n_x_n   };
    real_type * T0    { m_Tmat + iblock   * m_Tsize };
    integer   * P0    { m_Perm + iblock   * n       };
    real_type * work  { m_Fmat + n_thread * nr_x_nx };

    integer k{1};
    while ( k < nblk ) {

      real_type * Bjp { Bmat0             };
      real_type * Bj  { Bmat0 + k*n_x_nx  };
      real_type * Cjp { Cmat0             };
      real_type * Cj  { Cmat0 + k*nr_x_n  };
      real_type * Djp { Dmat0             };
      real_type * Dj  { Dmat0 + k*n_x_n   };
      real_type * Ejp { Emat0             };
      real_type * Ej  { Emat0 + k*n_x_n   };
      real_type * T   { T0    + k*m_Tsize };
      integer   * P   { P0    + k*n       };

      // -----------------------------------------------------------------------
      integer k_x_2{2*k};
      for ( integer j{iblock+k}; j < eblock; j += k_x_2 ) {

        GEcopy( n, n, Ejp, n, T,   n_x_2 ); // TOP
        GEcopy( n, n, Dj,  n, T+n, n_x_2 ); // BOTTOM

        if ( !buildT( n_thread, T, P ) ) {
          // factorization failed!
          m_ok_thread = false;
          return;
        }

        Zero_n( Dj,  n_x_n );
        Zero_n( Ejp, n_x_n );
        bool okk{applyT( n_thread, T, P, Djp, n, Dj, n, n ) &&
                 applyT( n_thread, T, P, Ejp, n, Ej, n, n )};
        if ( okk && m_nx > 0 ) okk = applyT( n_thread, T, P, Bjp, n, Bj, n, m_nx );

        if ( !okk ) {
          m_ok_thread = false;
          return;
        }

        if ( m_nr > 0 ) {

          if ( m_selected == BORDERED_Choice::QRP ) {
            integer i{n};
            do {
              --i;
              if ( P[i] > i ) alglin::swap( m_nr, Cj+i*m_nr, 1, Cj+P[i]*m_nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            SideMultiply::RIGHT,
            ULselect::UPPER,
            Transposition::NO,
            DiagonalType::NON_UNIT,
            m_nr, n, 1.0, T, n_x_2, Cj, m_nr
          );

          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Dj,  n,
             1.0, Cjp, m_nr
          );

          integer     jpp{ min(j+k,eblock) };
          real_type * Cpp{ m_Cmat + jpp*nr_x_n };

          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Ej,  n,
             1.0, Cpp, m_nr
          );
        }

        if ( nr_x_nx > 0 )
          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, m_nx, n,
            -1.0, Cj,   m_nr,
                  Bj,   n,
             1.0, work, m_nr // solo accumulo (F_mat)
          );

        // NEXT STEP
        T   += k_x_2*m_Tsize;
        P   += k_x_2*n;
        Djp += k_x_2*n_x_n;
        Dj  += k_x_2*n_x_n;
        Ejp += k_x_2*n_x_n;
        Ej  += k_x_2*n_x_n;
        Bj  += k_x_2*n_x_nx;
        Bjp += k_x_2*n_x_nx;
        Cj  += k_x_2*nr_x_n;
        Cjp += k_x_2*nr_x_n;
      }
      k *= 2;
    }
    m_kBlock[n_thread] = k;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::factorize_reduced() {

    integer const & n{m_block_size};
    integer k{1};
    while ( k < m_reduced_nblk ) {
      // -----------------------------------------------------------------------
      for ( integer jj{k}; jj < m_reduced_nblk; jj += 2*k ) {
        integer j  { m_iBlock[jj] };
        integer jp { m_iBlock[jj-k] };
        real_type * T   { m_Tmat + j*m_Tsize };
        integer   * P   { m_Perm + j*n       };
        real_type * Djp { m_Dmat + jp*n_x_n  };
        real_type * Dj  { m_Dmat + j*n_x_n   };
        real_type * Ejp { m_Emat + jp*n_x_n  };
        real_type * Ej  { m_Emat + j*n_x_n   };

        // check if factorization failed!
        GEcopy( n, n, Ejp, n, T,   n_x_2 ); // TOP
        GEcopy( n, n, Dj,  n, T+n, n_x_2 ); // BOTTOM
        if ( !buildT( 0, T, P ) ) return false;

        Zero_n( Dj, n_x_n );
        Zero_n( Ejp, n_x_n );
        bool okk{ applyT( 0, T, P, Djp, n, Dj, n, n ) &&
                  applyT( 0, T, P, Ejp, n, Ej, n, n ) };
        if ( okk && m_nx > 0 ) {
          real_type * Bj  { m_Bmat + j*n_x_nx  };
          real_type * Bjp { m_Bmat + jp*n_x_nx };
          okk = applyT( 0, T, P, Bjp, n, Bj, n, m_nx );
        }
        if ( !okk ) return false;

        if ( m_nr > 0 ) {
          real_type * Cj  { m_Cmat + j*nr_x_n  };
          real_type * Cjp { m_Cmat + jp*nr_x_n };
          if ( m_selected == BORDERED_Choice::QRP ) {
            integer i{n};
            do {
              --i;
              if ( P[i] > i ) alglin::swap( m_nr, Cj+i*m_nr, 1, Cj+P[i]*m_nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            SideMultiply::RIGHT,
            ULselect::UPPER,
            Transposition::NO,
            DiagonalType::NON_UNIT,
            m_nr, n, 1.0, T, n_x_2, Cj, m_nr
          );

          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Dj,  n,
             1.0, Cjp, m_nr
          );

          integer     jpp { m_iBlock[min(jj+k,m_reduced_nblk)] };
          real_type * Cpp { m_Cmat + jpp*nr_x_n };

          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Ej,  n,
             1.0, Cpp, m_nr
          );

        }

        if ( nr_x_nx > 0 ) {
          real_type * Cj{ m_Cmat + j*nr_x_n };
          real_type * Bj{ m_Bmat + j*n_x_nx };
          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, m_nx, n,
            -1.0, Cj,     m_nr,
                  Bj,     n,
             1.0, m_Fmat, m_nr
          );
        }
      }
      k *= 2;
    }
    return true;
  }

  /*\
   |   _                 _                   _
   |  | | ___   __ _  __| |   __ _ _ __   __| |
   |  | |/ _ \ / _` |/ _` |  / _` | '_ \ / _` |
   |  | | (_) | (_| | (_| | | (_| | | | | (_| |
   |  |_|\___/ \__,_|\__,_|  \__,_|_| |_|\__,_|
   |    __            _             _           _           _
   |   / _| __ _  ___| |_ ___  _ __(_)_______  | | __ _ ___| |_
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \ | |/ _` / __| __|
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/ | | (_| \__ \ |_
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___| |_|\__,_|___/\__|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::load_and_factorize_last() {
    integer const & nblock { m_number_of_blocks };
    integer const & n      { m_block_size };
   /*
    //    n   n  qx  nx
    //  +---+---+---+---+
    //  | D | E | 0 | B | n
    //  +===+===+===+===+
    //  |   |   |   |   |
    //  |H0 |HN |Hq |Hp | n+qr
    //  |   |   |   |   |
    //  +---+---+---+---+
    //  | C | C |Cq | F | nr
    //  +---+---+---+---+
    */
    real_type * Cnb { m_Cmat + nblock*nr_x_n };
    real_type * W0  { m_Hmat };
    real_type * WN  { W0+n*m_Nr };
    real_type * Wq  { WN+n*m_Nr };
    real_type * Wp  { Wq+m_qx*m_Nr };

    GEcopy( n, n,    m_Dmat, n, W0, m_Nr );
    GEcopy( n, n,    m_Emat, n, WN, m_Nr );
    GEzero( n, m_qx,            Wq, m_Nr );
    GEcopy( n, m_nx, m_Bmat, n, Wp, m_Nr );

    GEcopy( n+m_qr, m_Nc, m_H0Nqp, n+m_qr, m_Hmat+n, m_Nr );

    integer offs { n_x_2+m_qr };

    GEcopy( m_nr, n,    m_Cmat,  m_nr, W0+offs, m_Nr );
    GEcopy( m_nr, n,    Cnb,     m_nr, WN+offs, m_Nr );
    GEcopy( m_nr, m_qx, m_Cqmat, m_nr, Wq+offs, m_Nr );
    GEcopy( m_nr, m_nx, m_Fmat,  m_nr, Wp+offs, m_Nr );

    bool ok = false;
    switch ( m_last_block_selected ) {
    case BORDERED_LAST_Choice::LU:
      ok = m_last_lu.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::LUPQ:
      ok = m_last_lupq.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::QR:
      ok = m_last_qr.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::QRP:
      ok = m_last_qrp.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::SVD:
      ok = m_last_svd.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::LSS:
      ok = m_last_lss.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::LSY:
      ok = m_last_lsy.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_Choice::PINV:
      // try to factorize using PINV
      if ( m_Nc > m_Nr ) ok = m_last_pinv.t_factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      else               ok = m_last_pinv.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    }
    if ( ok ) {
      m_last_block_use_solver2 = false;
    } else {
      m_last_block_use_solver2 = true;
      switch ( m_last_block_selected2 ) {
      case BORDERED_LAST_Choice2::NONE:
        break;
      case BORDERED_LAST_Choice2::SVD:
        ok = m_last_svd.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
        break;
      case BORDERED_LAST_Choice2::LSS:
        ok = m_last_lss.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
        break;
      case BORDERED_LAST_Choice2::LSY:
        ok = m_last_lsy.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
        break;
      case BORDERED_LAST_Choice2::PINV:
        // try to factorize using PINV
        if ( m_Nc > m_Nr ) ok = m_last_pinv.t_factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
        else               ok = m_last_pinv.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
        break;
      }
    }
    if ( !ok )
      m_last_error = "BorderedCR<t_Value>::load_and_factorize_last failed\n";

    return ok;
  }

  /*\
   |             _               _           _
   |   ___  ___ | |_   _____    | | __ _ ___| |_
   |  / __|/ _ \| \ \ / / _ \   | |/ _` / __| __|
   |  \__ \ (_) | |\ V /  __/   | | (_| \__ \ |_
   |  |___/\___/|_| \_/ \___|___|_|\__,_|___/\__|
   |                       |_____|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_last( real_type x[] ) const {

    integer const & nblock { m_number_of_blocks };
    integer const & n      { m_block_size };
    real_type * X{ x + (nblock-1)*n };
    // sposto primo blocco rhs in fondo
    swap( n, X, 1, x, 1 ); // uso x stesso come temporaneo
    bool ok{false};
    if ( m_last_block_use_solver2 ) {
      switch ( m_last_block_selected2 ) {
        case BORDERED_LAST_Choice2::NONE: break;
        case BORDERED_LAST_Choice2::SVD:  ok = m_last_svd.solve( X );  break;
        case BORDERED_LAST_Choice2::LSS:  ok = m_last_lss.solve( X );  break;
        case BORDERED_LAST_Choice2::LSY:  ok = m_last_lsy.solve( X );  break;
        case BORDERED_LAST_Choice2::PINV:
          if ( m_Nc > m_Nr ) ok = m_last_pinv.t_mult_inv( X, 1, X, 1 );
          else               ok = m_last_pinv.mult_inv( X, 1, X, 1 );
          break;
      }
    } else {
      switch ( m_last_block_selected ) {
        case BORDERED_LAST_Choice::LU:   ok = m_last_lu.solve( X );   break;
        case BORDERED_LAST_Choice::LUPQ: ok = m_last_lupq.solve( X ); break;
        case BORDERED_LAST_Choice::QR:   ok = m_last_qr.solve( X );   break;
        case BORDERED_LAST_Choice::QRP:  ok = m_last_qrp.solve( X );  break;
        case BORDERED_LAST_Choice::SVD:  ok = m_last_svd.solve( X );  break;
        case BORDERED_LAST_Choice::LSS:  ok = m_last_lss.solve( X );  break;
        case BORDERED_LAST_Choice::LSY:  ok = m_last_lsy.solve( X );  break;
        case BORDERED_LAST_Choice::PINV:
          if ( m_Nc > m_Nr ) ok = m_last_pinv.t_mult_inv( X, 1, X, 1 );
          else               ok = m_last_pinv.mult_inv( X, 1, X, 1 );
          break;
      }
    }
    if ( ok ) alglin::swap( n, X, 1, x, 1 );
    else      m_last_error = "BorderedCR<t_Value>::solve_last( x ) failed\n";
    return ok;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_last(
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {

    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};
    real_type * X{x + (nblock-1)*n};
    // sposto primo blocco rhs in fondo
    for ( integer i{0}; i < nrhs; ++i )
      alglin::swap( n, X+i*ldX, 1, x+i*ldX, 1 );
    bool ok = false;
    if ( m_last_block_use_solver2 ) {
      switch ( m_last_block_selected2 ) {
        case BORDERED_LAST_Choice2::NONE: break;
        case BORDERED_LAST_Choice2::SVD:  ok = m_last_svd.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice2::LSS:  ok = m_last_lss.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice2::LSY:  ok = m_last_lsy.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice2::PINV:
          if ( m_Nc > m_Nr ) ok = m_last_pinv.t_mult_inv( nrhs, X, ldX, X, ldX );
          else               ok = m_last_pinv.mult_inv( nrhs, X, ldX, X, ldX );
          break;
      }
    } else {
      switch ( m_last_block_selected ) {
        case BORDERED_LAST_Choice::LU:   ok = m_last_lu.solve( nrhs, X, ldX  );   break;
        case BORDERED_LAST_Choice::LUPQ: ok = m_last_lupq.solve( nrhs, X, ldX  ); break;
        case BORDERED_LAST_Choice::QR:   ok = m_last_qr.solve( nrhs, X, ldX );   break;
        case BORDERED_LAST_Choice::QRP:  ok = m_last_qrp.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice::SVD:  ok = m_last_svd.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice::LSS:  ok = m_last_lss.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice::LSY:  ok = m_last_lsy.solve( nrhs, X, ldX  );  break;
        case BORDERED_LAST_Choice::PINV:
          if ( m_Nc > m_Nr ) ok = m_last_pinv.t_mult_inv( nrhs, X, ldX, X, ldX );
          else               ok = m_last_pinv.mult_inv( nrhs, X, ldX, X, ldX );
          break;
      }
    }
    if ( ok )
      for ( integer i{0}; i < nrhs; ++i )
        alglin::swap( n, X+i*ldX, 1, x+i*ldX, 1 );
    if ( !ok )
      m_last_error = "BorderedCR<t_Value>::solve_last( nrhs, x, ldX ) failed\n";
    return ok;
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve( real_type x[] ) const {
    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    integer const & npb    {m_used_parallel_block};
    real_type * xb{x + (nblock+1)*n + m_qr}; // deve essere b!
    if ( npb > 1 ) {
      if ( m_solve_use_thread ) {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          m_TP->run( &BorderedCR<t_Value>::forward, this, n_thread, x );
        m_TP->wait();
      } else {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          forward( n_thread, x );
      }
      if ( m_nr > 0 ) {
        MapVector<t_Value> XB( xb, m_nr );
        for ( integer n_thread{0}; n_thread < npb; ++n_thread ) {
          MapVector<t_Value> WORK( m_xb_thread[size_t(n_thread)], m_nr );
          XB.noalias() += WORK;
        }
      }
      if ( !forward_reduced( x, xb ) ) return false;
    } else {
      forward( 0, x );
      if ( m_nr > 0 ) {
        MapVector<t_Value> XB( xb, m_nr );
        MapVector<t_Value> WORK( m_xb_thread[0], m_nr );
        XB.noalias() += WORK;
      }
    }

    if ( !solve_last( x ) ) return false;

    if ( npb > 1 ) {
      backward_reduced(x);
      if ( m_solve_use_thread ) {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          m_TP->run( &BorderedCR<t_Value>::backward, this, n_thread, x );
        m_TP->wait();
      } else {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          backward( n_thread, x );
      }
    } else {
      backward( 0, x );
    }
    return true;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve(
    integer   nrhs,
    real_type rhs[],
    integer   ldRhs
  ) const {
    integer const & npb    {m_used_parallel_block};
    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    MapMatrix<t_Value> XB( rhs + (nblock+1)*n + m_qr, ldRhs, nrhs );

    m_ok_thread = true;

    integer nn{m_nr*nrhs};
    real_type * work{m_work_mem.realloc( npb * nn )};

    if ( npb > 1 ) {
      if ( m_solve_use_thread ) {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          m_TP->run( &BorderedCR<t_Value>::forward_n, this, n_thread, nrhs, rhs, ldRhs, work + n_thread * nn );
        m_TP->wait();
      } else {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          forward_n( n_thread, nrhs, rhs, ldRhs, work + n_thread * nn );
      }
      if ( m_nr > 0 ) {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread ) {
          MapMatrix<t_Value> WORK( work + n_thread * nn, m_nr, nrhs );
          XB.topRows( m_nr ).noalias() += WORK;
        }
      }
      if ( m_ok_thread ) m_ok_thread = forward_n_reduced( nrhs, rhs, ldRhs );
    } else {
      forward_n( 0, nrhs, rhs, ldRhs, work );
      if ( m_nr > 0 ) {
        MapMatrix<t_Value> WORK( work, m_nr, nrhs );
        XB.topRows( m_nr ).noalias() += WORK;
      }
    }

    if ( m_ok_thread ) m_ok_thread = solve_last( nrhs, rhs, ldRhs );
    if ( !m_ok_thread ) return false;

    if ( npb > 1 ) {
      backward_n_reduced( nrhs, rhs, ldRhs );
      if ( m_solve_use_thread ) {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          m_TP->run( &BorderedCR<t_Value>::backward_n, this, n_thread, nrhs, rhs, ldRhs );
        m_TP->wait();
      } else {
        for ( integer n_thread{0}; n_thread < npb; ++n_thread )
          backward_n( n_thread, nrhs, rhs, ldRhs );
      }
    } else {
      backward_n( 0, nrhs, rhs, ldRhs );
    }
    return true;
  }

  /*\
   |    __                                  _
   |   / _| ___  _ ____      ____ _ _ __ __| |
   |  | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |
   |  |  _| (_) | |   \ V  V / (_| | | | (_| |
   |  |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::forward( integer n_thread, real_type x[] ) const {
    integer const & n = m_block_size;

    integer   iblock { m_iBlock[2*n_thread+0]        };
    integer   eblock { m_iBlock[2*n_thread+1]        };
    integer   nblk   { eblock - iblock               };
    real_type * x0   { x      + iblock*n             };
    real_type * T0   { m_Tmat + iblock*m_Tsize       };
    integer   * P0   { m_Perm + iblock*n             };
    real_type * C0   { m_Cmat + iblock*nr_x_n        };
    real_type * work { m_xb_thread[size_t(n_thread)] };

    // se serve accumula in work
    if ( m_nr > 0 ) Zero_n( work, m_nr );

    integer k = 1;
    while ( k < nblk ) {
      real_type * xj  { x0 + k*n        };
      real_type * xjp { x0              };
      real_type * T   { T0 + k*m_Tsize  };
      integer   * P   { P0 + k*n        };
      real_type * Cj  { C0 + k*nr_x_n   };
      integer   k_x_2 { 2*k             };
      for ( integer jj{k}; jj < nblk; jj += k_x_2 ) {
        if ( !applyT( n_thread, T, P, xjp, xj ) ) return false;
        if ( m_nr > 0 )
          gemv(
            Transposition::NO,
            m_nr, n, -1.0,
            Cj, m_nr, xj, 1,
            1.0, work, 1
          ); // solo accumulato
        xj  += k_x_2*n;
        xjp += k_x_2*n;
        T   += k_x_2*m_Tsize;
        P   += k_x_2*n;
        Cj  += k_x_2*nr_x_n;
      }
      k *= 2;
    }

    return true;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_n(
    integer   n_thread,
    integer   nrhs,
    real_type x[],
    integer   ldX,
    real_type work[]
  ) const {
    integer const & n{m_block_size};
    integer iblock { m_iBlock[2*n_thread+0]  };
    integer eblock { m_iBlock[2*n_thread+1]  };
    integer nblk   { eblock - iblock         };
    real_type * x0 { x      + iblock*n       };
    real_type * T0 { m_Tmat + iblock*m_Tsize };
    integer   * P0 { m_Perm + iblock*n       };
    real_type * C0 { m_Cmat + iblock*nr_x_n  };

    // se serve accumula in work
    if ( m_nr > 0 ) Zero_n( work, m_nr*nrhs );

    integer k{1};
    while ( k < nblk ) {
      real_type * xj  { x0 + k*n       };
      real_type * xjp { x0             };
      real_type * T   { T0 + k*m_Tsize };
      integer   * P   { P0 + k*n       };
      real_type * Cj  { C0 + k*nr_x_n  };
      integer   k_x_2 { 2*k            };

      for ( integer jj{k}; jj < nblk; jj += k_x_2 ) {
        if ( !applyT( n_thread, T, P, xjp, ldX, xj, ldX, nrhs ) ) {
          m_ok_thread = false;
          return;
        }
        if ( m_nr > 0 ) {
          // accumula in work
          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, nrhs, n,
            -1.0, Cj, m_nr,
                  xj, ldX,
             1.0, work, m_nr
          );
        }
        xj  += k_x_2*n;
        xjp += k_x_2*n;
        T   += k_x_2*m_Tsize;
        P   += k_x_2*n;
        Cj  += k_x_2*nr_x_n;
      }
      k *= 2;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::forward_reduced(
    real_type x[],
    real_type xb[]
  ) const {
    integer const & n{m_block_size};

    integer k{1};
    while ( k < m_reduced_nblk ) {
      for ( integer jj = k; jj < m_reduced_nblk; jj += 2*k ) {
        integer j{m_iBlock[jj]};
        integer jp{m_iBlock[jj-k]};
        real_type const * T   { m_Tmat + j*m_Tsize };
        integer   const * P   { m_Perm + j*n       };
        real_type       * xj  { x + j*n            };
        real_type       * xjp { x + jp*n           };
        if ( !applyT( 0, T, P, xjp, xj ) ) return false;
        if ( m_nr > 0 ) {
          real_type * Cj = m_Cmat + j*nr_x_n;
          alglin::gemv(
            Transposition::NO, m_nr, n,
            -1.0, Cj, m_nr, xj, 1, 1.0, xb, 1
          );
        }
      }
      k *= 2;
    }
    return true;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::forward_n_reduced(
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    real_type * xb{ x + (nblock+1)*n + m_qr };
    integer k{1};
    while ( k < m_reduced_nblk ) {
      for ( integer jj{k}; jj < m_reduced_nblk; jj += 2*k ) {
        integer j{m_iBlock[jj]};
        integer jp{m_iBlock[jj-k]};
        real_type const * T   { m_Tmat + j*m_Tsize };
        integer   const * P   { m_Perm + j*n       };
        real_type       * xj  { x + j*n            };
        real_type       * xjp { x + jp*n           };
        if ( !applyT( 0, T, P, xjp, ldX, xj, ldX, nrhs ) ) return false;
        if ( m_nr > 0 ) {
          // per ora non serve
          // #pragma omp critical
          real_type * Cj{m_Cmat + j*nr_x_n};
          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            m_nr, nrhs, n,
            -1.0, Cj, m_nr,
                  xj, ldX,
             1.0, xb, ldX
          );
        }
      }
      k *= 2;
    }
    return true;
  }

  /*\
   |   _                _                           _
   |  | |__   __ _  ___| | ____      ____ _ _ __ __| |
   |  | '_ \ / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
   |  | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
   |  |_.__/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward( integer n_thread, real_type x[] ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    real_type * xn { x + (nblock+1)*n + m_qx };
    integer iblock { m_iBlock[2*n_thread+0]  };
    integer eblock { m_iBlock[2*n_thread+1]  };
    real_type * x0 { x      + iblock*n       };
    real_type * B0 { m_Bmat + iblock*n_x_nx  };
    real_type * D0 { m_Dmat + iblock*n_x_n   };
    real_type * E0 { m_Emat + iblock*n_x_n   };
    real_type * T0 { m_Tmat + iblock*m_Tsize };
    integer k = m_kBlock[n_thread];
    while ( (k/=2) > 0 ) {
      real_type * xj { x0 + k*n       };
      real_type * xp { x0             };
      real_type * Bj { B0 + k*n_x_nx  };
      real_type * Dj { D0 + k*n_x_n   };
      real_type * Ej { E0 + k*n_x_n   };
      real_type * T  { T0 + k*m_Tsize };
      integer k_x_2  { 2*k            };
      for ( integer j{iblock+k}; j < eblock; j += k_x_2 ) {
        integer     jpp{min(j+k,eblock)};
        real_type * xpp{x + jpp*n};
        alglin::gemv(
          Transposition::NO,
          n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1
        );
        alglin::gemv(
          Transposition::NO,
          n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1
        );
        if ( m_nx > 0 )
          gemv(
            Transposition::NO,
            n, m_nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1
          );
        alglin::trsv(
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          n, T, n_x_2, xj, 1
        );
        if ( m_selected == BORDERED_Choice::QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              alglin::swap( xj[i], xj[P[i]] );
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx;
        Dj += k_x_2*n_x_n;
        Ej += k_x_2*n_x_n;
        xj += k_x_2*n;
        xp += k_x_2*n;
        T  += k_x_2*m_Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n(
    integer   n_thread,
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {

    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    real_type * xn { x + (nblock+1)*n + m_qx };
    integer iblock { m_iBlock[2*n_thread+0]  };
    integer eblock { m_iBlock[2*n_thread+1]  };
    real_type * x0 { x      + iblock*n       };
    real_type * B0 { m_Bmat + iblock*n_x_nx  };
    real_type * D0 { m_Dmat + iblock*n_x_n   };
    real_type * E0 { m_Emat + iblock*n_x_n   };
    real_type * T0 { m_Tmat + iblock*m_Tsize };
    integer k{m_kBlock[n_thread]};
    while ( (k/=2) > 0 ) {
      real_type * xj { x0 + k*n       };
      real_type * xp { x0             };
      real_type * Bj { B0 + k*n_x_nx  };
      real_type * Dj { D0 + k*n_x_n   };
      real_type * Ej { E0 + k*n_x_n   };
      real_type * T  { T0 + k*m_Tsize };
      integer  k_x_2 { 2*k            };
      for ( integer j{iblock+k}; j < eblock; j += k_x_2 ) {
        integer     jpp{min(j+k,eblock)};
        real_type * xpp{x + jpp*n};
        alglin::gemm(
          Transposition::NO,
          Transposition::NO,
          n, nrhs, n,
          -1.0, Dj,  n,
                xp,  ldX,
           1.0, xj,  ldX
        );
        alglin::gemm(
          Transposition::NO,
          Transposition::NO,
          n, nrhs, n,
          -1.0, Ej,  n,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( m_nx > 0 )
          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            n, nrhs, m_nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        alglin::trsm(
          SideMultiply::LEFT,
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          n, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( m_selected == BORDERED_Choice::QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              alglin::swap( nrhs, xj+i, ldX, xj+P[i], ldX );
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx;
        Dj += k_x_2*n_x_n;
        Ej += k_x_2*n_x_n;
        xj += k_x_2*n;
        xp += k_x_2*n;
        T  += k_x_2*m_Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_reduced( real_type x[] ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    real_type * xn{ x + (nblock+1)*n + m_qx };
    integer k{1};
    while ( k < m_reduced_nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj{k}; jj < m_reduced_nblk; jj += 2*k ) {
        integer     j   { m_iBlock[jj] };
        integer     jp  { m_iBlock[jj-k] };
        integer     jpp { m_iBlock[min(jj+k,m_reduced_nblk)] };
        real_type * Dj  { m_Dmat + j*n_x_n };
        real_type * Ej  { m_Emat + j*n_x_n };
        real_type * xj  { x + j*n };
        real_type * xp  { x + jp*n };
        real_type * xpp { x + jpp*n };
        alglin::gemv(
          Transposition::NO,
          n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1
        );
        alglin::gemv(
          Transposition::NO,
          n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1
        );
        if ( m_nx > 0 ) {
          real_type * Bj = m_Bmat + j*n_x_nx;
          gemv(
            Transposition::NO,
            n, m_nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1
          );
        }
        real_type const * T{ m_Tmat + j*m_Tsize };
        alglin::trsv(
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          n, T, n_x_2, xj, 1
        );
        if ( m_selected == BORDERED_Choice::QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] );
        }
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n_reduced(
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {

    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    real_type * xn{x + (nblock+1)*n + m_qx};
    integer k{1};
    while ( k < m_reduced_nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj{k}; jj < m_reduced_nblk; jj += 2*k ) {
        integer     j   { m_iBlock[jj] };
        integer     jp  { m_iBlock[jj-k] };
        integer     jpp { m_iBlock[min(jj+k,m_reduced_nblk)] };
        real_type * Dj  { m_Dmat + j*n_x_n };
        real_type * Ej  { m_Emat + j*n_x_n };
        real_type * xj  { x + j*n };
        real_type * xjp { x + jp*n };
        real_type * xpp { x + jpp*n };
        alglin::gemm(
          Transposition::NO,
          Transposition::NO,
          n, nrhs, n,
          -1.0, Dj,  n,
                xjp, ldX,
           1.0, xj,  ldX
        );
        alglin::gemm(
          Transposition::NO,
          Transposition::NO,
          n, nrhs, n,
          -1.0, Ej,  n,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( m_nx > 0 ) {
          real_type * Bj = m_Bmat + j*n_x_nx;
          alglin::gemm(
            Transposition::NO,
            Transposition::NO,
            n, nrhs, m_nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        real_type const * T{ m_Tmat + j*m_Tsize };
        trsm(
          SideMultiply::LEFT,
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          n, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( m_selected == BORDERED_Choice::QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              alglin::swap( nrhs, xj+i, ldX, xj+P[i], ldX );
        }
      }
    }
  }

  /*\
   |  _  _  ___ _____  __   _____ _____
   |  | \| |/ _ \_   _| \ \ / / __|_   _|
   |  | .` | (_) || |    \ V /| _|  | |
   |  |_|\_|\___/ |_|     |_| |___| |_|
   |
   |   ___ __  __ ___ _    ___ __  __ ___ _  _ _____ ___ ___
   |  |_ _|  \/  | _ \ |  | __|  \/  | __| \| |_   _| __|   \
   |   | || |\/| |  _/ |__| _|| |\/| | _|| .` | | | | _|| |) |
   |  |___|_|  |_|_| |____|___|_|  |_|___|_|\_| |_| |___|___/
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::t_solve( real_type [] ) const {
    UTILS_ERROR0( "BorderedCR::t_solve() not defined\n" );
    return false;
  }

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::t_solve( integer, real_type [], integer ) const {
    UTILS_ERROR0( "BorderedCR::t_solve() not defined\n" );
    return false;
  }

  /*\
   |             _     _ __  __
   |    __ _  __| | __| |  \/  |_   __
   |   / _` |/ _` |/ _` | |\/| \ \ / /
   |  | (_| | (_| | (_| | |  | |\ V /
   |   \__,_|\__,_|\__,_|_|  |_| \_/
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::Mv( real_type const x[], real_type res[] ) const {
    Zero_n( res, nrows() );
    add_Mv( 1.0, x, res );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::add_Mv(
    real_type       alpha,
    real_type const x[],
    real_type       res[]
  ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    // internal blocks block
    t_Value const * D  { m_Dmat        };
    t_Value const * E  { m_Emat        };
    t_Value const * B  { m_Bmat        };
    t_Value const * xx { x             };
    t_Value const * xe { x  + nblock*n };
    t_Value const * xq { xe + n        };
    t_Value const * xb { xq + m_qx     };
    t_Value *       yy = res;
    for ( integer i{0}; i < nblock; ++i ) {
      alglin::gemv(
        Transposition::NO, n, n, alpha, D, n, xx, 1, 1.0, yy, 1
      );
      xx += n;
      alglin::gemv(
        Transposition::NO, n, n, alpha, E, n, xx, 1, 1.0, yy, 1
      );
      alglin::gemv(
        Transposition::NO, n, m_nx, alpha, B, n, xb, 1, 1.0, yy, 1
      );
      yy += n;
      D  += n_x_n;
      E  += n_x_n;
      B  += n_x_nx;
    }

    integer     m{n+m_qr};
    real_type * H{m_H0Nqp};
    alglin::gemv( Transposition::NO, m, n,    alpha, H, m, x,  1, 1.0, yy, 1 ); H += m * n;
    alglin::gemv( Transposition::NO, m, n,    alpha, H, m, xe, 1, 1.0, yy, 1 ); H += m * n;
    alglin::gemv( Transposition::NO, m, m_qx, alpha, H, m, xq, 1, 1.0, yy, 1 ); H += m * m_qx;
    alglin::gemv( Transposition::NO, m, m_nx, alpha, H, m, xb, 1, 1.0, yy, 1 );

    if ( m_nr > 0 ) {
      yy += m;
      alglin::gemv(
        Transposition::NO, m_nr, m_nx,
        alpha, m_Fmat, m_nr, xb, 1, 1.0, yy, 1
      );
      t_Value const * C = m_Cmat;
      xx = x;
      for ( integer i = 0; i <= nblock; ++i ) {
        alglin::gemv(
          Transposition::NO, m_nr, n,
          alpha, C, m_nr, xx, 1, 1.0, yy, 1
        );
        xx += n; C += nr_x_n;
      }
      alglin::gemv(
        Transposition::NO, m_nr, m_qx,
        alpha, m_Cqmat, m_nr, xx, 1, 1.0, yy, 1
      );
    }
  }

  /*\
   |  ____
   | / ___| _ __   __ _ _ __ ___  ___
   | \___ \| '_ \ / _` | '__/ __|/ _ \
   |  ___) | |_) | (_| | |  \__ \  __/
   | |____/| .__/ \__,_|_|  |___/\___|
   |       |_|
   |
  \*/

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_B(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    integer i0{nbl*n + offs};
    integer j0{(nblock+1)*n + m_qx + offs};
    for ( integer ij{0}; ij < n_x_nx; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_B( integer nbl, real_type V[] ) const {
    Copy_n( m_Bmat + nbl*n_x_nx, n_x_nx, V );
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_C(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    integer i0{(nblock+1)*n + m_qr + offs};
    integer j0{nbl*n + offs};
    for ( integer ij{0}; ij < nr_x_n; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_C( integer nbl, real_type V[] ) const {
    Copy_n( m_Cmat + nbl*nr_x_n, nr_x_n, V );
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_D(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & n{m_block_size};
    integer i0{nbl*n + offs};
    integer j0{nbl*n + offs};
    for ( integer ij{0}; ij < n_x_n; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_D( integer nbl, real_type V[] ) const {
    Copy_n( m_Dmat + nbl*n_x_n, n_x_n, V );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_E(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & n{m_block_size};
    integer i0{nbl*n + offs};
    integer j0{(nbl+1)*n + offs};
    for ( integer ij{0}; ij < n_x_n; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_E( integer nbl, real_type V[] ) const {
    Copy_n( m_Emat + nbl*n_x_n, n_x_n, V );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_F(
    integer I[], integer J[], integer offs
  ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};
    integer i0{(nblock+1)*n + m_qr + offs};
    integer j0{(nblock+1)*n + m_qx + offs};
    for ( integer ij{0}; ij < nr_x_nx; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_F( real_type V[] ) const {
    Copy_n( m_Fmat, nr_x_nx, V );
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_Cq(
    integer I[], integer J[], integer offs
  ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    integer i0{(nblock+1)*n + m_qr + offs};
    integer j0{(nblock+1)*n + offs};
    for ( integer ij{0}; ij < nr_x_qx; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_Cq( real_type V[] ) const {
    Copy_n( m_Cqmat, nr_x_qx, V );
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::pattern_H( integer I[], integer J[], integer offs ) const {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};
    integer nqr{n + m_qr};
    integer nnz{nqr * m_Nc};
    integer i0{nblock*n};
    integer j0{i0 - n};
    for ( integer ij{0}; ij < nnz; ++ij ) {
      I[ij] = i0 + (ij % nqr) + offs;
      integer j = integer(ij/nqr); if ( j >= n ) j += j0;
      J[ij] = j + offs;
    }
    return nnz;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::values_H( real_type V[] ) const {
    integer const & n{m_block_size};
    integer nnz{(n + m_qr) * m_Nc};
    Copy_n( m_H0Nqp, nnz, V );
    return nnz;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparse_pattern(
    integer I[],
    integer J[],
    integer offs
  ) const {
    integer const & nblock{m_number_of_blocks};

    integer kkk{0};
    for ( integer nbl{0}; nbl < nblock; ++nbl ) {
      kkk += this->pattern_D( nbl, I+kkk, J+kkk, offs );
      kkk += this->pattern_E( nbl, I+kkk, J+kkk, offs );
      kkk += this->pattern_B( nbl, I+kkk, J+kkk, offs );
    }

    // H
    kkk += this->pattern_H( I+kkk, J+kkk, offs );

    // C
    for ( integer nbl = 0; nbl <= nblock; ++nbl )
      kkk += this->pattern_C( nbl, I+kkk, J+kkk, offs );

    // F
    kkk += this->pattern_F( I+kkk, J+kkk, offs );

    // Cq
    kkk += this->pattern_Cq( I+kkk, J+kkk, offs );

    if ( kkk == sparse_nnz() ) {
      BABD_LAST_ERROR(
        "BorderedCR::sparse_pattern( I, J, offs ), inserted {} values, expected {}\n",
        kkk, sparse_nnz()
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparse_values( real_type V[] ) const {
    integer const & nblock{m_number_of_blocks};
    integer kkk{0};
    for ( integer nbl{0}; nbl < nblock; ++nbl ) {
      kkk += this->values_D( nbl, V+kkk );
      kkk += this->values_E( nbl, V+kkk );
      kkk += this->values_B( nbl, V+kkk );
    }

    // H
    kkk += this->values_H( V+kkk );

    // C
    for ( integer nbl{0}; nbl <= nblock; ++nbl )
      kkk += this->values_C( nbl, V+kkk );

    // F
    kkk += this->values_F( V+kkk );

    // Cq
    kkk += this->values_Cq( V+kkk );

    if ( kkk == sparse_nnz() ) {
      BABD_LAST_ERROR(
        "BorderedCR::sparse_values( V ), inserted {} values, expected {}\n",
        kkk, sparse_nnz()
      );
      Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparse_load(
    real_type const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz
  ) {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};

    integer const rH   {n*nblock};
    integer const rC   {rH + n + m_qr};
    integer const nrow {rC + m_nr};

    integer const cCq  {rH  + n};
    integer const cF   {cCq + m_qx};
    integer const ncol {cF  + m_nx};

    fill_zero();

    for ( integer kkk{0}; kkk < M_nnz; ++kkk ) {
      integer   i{M_row[kkk] - r_offs};
      integer   j{M_col[kkk] - c_offs};
      real_type v{M_values[kkk]};
      // cerca blocco
      bool ok{true};
      if ( i < rH ) {
        if ( j < cCq ) { // DE
          // cerca blocchi
          integer ib = i/n;
          integer jb = j/n;
          if ( ib == jb ) {
            D(ib,i%n,j%n) = v;
          } else if ( ib+1 == jb ) {
            E(ib,i%n,j%n) = v;
          } else {
            ok = false;
          }
        } else if ( j < cF ) { // Hq
          ok = false;
        } else if ( j < ncol ) { // B
          integer ib = i/n;
          B(ib,i%n,j-cF) = v;
        } else {
          ok = false;
        }
      } else if ( i < rC ) {
        if ( j < n ) { // H0
          H(i-rH,j) = v;
        } else if ( j < rH ) {
          ok = false;
        } else if ( j < ncol ) { // HN,Hq,Hp
          H(i-rH,j-rH+n) = v;
        } else {
          ok = false;
        }
      } else if ( i < nrow ) {
        if ( j < cCq ) {
          integer jb = j/n;
          C(jb,i-rC,j%n) = v;
        } else if ( j < cF ) {
          Cq(i-rC,j-cCq) =v;
        } else if ( j < ncol ) {
          F(i-rC,j-cF) = v;
        } else {
          ok = false;
        }
      } else {
        ok = false;
      }
      if ( !ok ) {
        BABD_LAST_ERROR(
          "in BorderedCR<t_Value>::sparse_load, "
          "indices (i,j) = ( {}, {}) out of pattern!\n",
          M_row[kkk], M_col[kkk]
        );
        Utils::Runtime_Error( m_last_error, __FILE__, __LINE__ );
      }
    }
  }

  // ---------------------------------------------------------------------------

  /*\
   |   __  __   _ _____ _      _   ___
   |  |  \/  | /_\_   _| |    /_\ | _ )
   |  | |\/| |/ _ \| | | |__ / _ \| _ \
   |  |_|  |_/_/ \_\_| |____/_/ \_\___/
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::print_matlab_script( ostream_type & stream ) const {
    integer nnz = this->sparse_nnz();
    Malloc<t_Value> mem(" BorderedCR::print_matlab_script real");
    Malloc<integer> memi(" BorderedCR::print_matlab_script integer");
    memi.allocate( size_t(2*nnz) );
    t_Value * V{mem.malloc( size_t(nnz) )};
    integer * I{memi( size_t(nnz) )};
    integer * J{memi( size_t(nnz) )};

    this->sparse_pattern( I, J, 1 );
    this->sparse_values( V );

    fmt::print( stream, "I = [ {}", I[0] );
    for ( integer i{1}; i < nnz; ++i ) {
      if ( (i % 20) == 0 ) fmt::print( stream, ", ...\n  {}", I[i] );
      else                 fmt::print( stream, ",  {}", I[i] );
    }
    fmt::print( stream, " ];\n\n" );

    fmt::print( stream, "J = [ {}", J[0] );
    for ( integer i{1}; i < nnz; ++i ) {
      if ( (i % 20) == 0 ) fmt::print( stream, ", ...\n  {}", J[i] );
      else                 fmt::print( stream, ",  {}", J[i] );
    }
    fmt::print( stream, " ];\n\n" );

    fmt::print( stream, "V = [ {}", V[0] );
    for ( integer i{1}; i < nnz; ++i ) {
      if ( (i % 20) == 0 ) fmt::print( stream, ", ...\n  {}", V[i] );
      else                 fmt::print( stream, ",  {}", V[i] );
    }
    fmt::print( stream,
      " ];\n\n"
      "nr  = {};\n"
      "nc  = {};\n\n"
      "MAT = sparse( I, J, V, nr, nc );\n",
      this->nrows(), this->ncols()
    );
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::check_matrix() {

    try{

      integer const & nblock{m_number_of_blocks};
      integer const & n{m_block_size};
      {
        integer NR{n+m_qr};
        integer NC{m_Nc};
        integer NN{NR*NC};
        for( integer i{0}; i < NN; ++i ) {
          integer r{ i % NR };
          integer c{ integer(i / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_H0Nqp[i]),
            "BorderedCR::check_matrix, found {} at position ({},{}) of {}x{} matrix H0Hqp\n",
            m_H0Nqp[i], r, c, NR, NC
         );
        }
      }
      {
        integer NR{n};
        integer NC{m_nx};
        integer SB{NR*NC};
        integer NN{SB*nblock};
        for( integer i{0}; i < NN; ++i ) {
          integer nb{ integer(i/SB) };
          integer r{ i % NR };
          integer c{ integer((i-nb*SB) / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_Bmat[i]),
            "BorderedCR::check_matrix, found {} at block {} position ({},{}) of {}x{} matrix B\n",
            m_Bmat[i], nb, r, c, NR, NC
         );
        }
      }
      {
        integer NR{m_nr};
        integer NC{n};
        integer SB{NR*NC};
        integer NN{SB*(nblock+1)};
        for( integer i{0}; i < NN; ++i ) {
          integer nb{ integer(i/SB) };
          integer r{ i % NR };
          integer c{ integer((i-nb*SB) / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_Cmat[i]),
            "BorderedCR::check_matrix, found {} at block {} position ({},{}) of {}x{} matrix C\n",
            m_Cmat[i], nb, r, c, NR, NC
         );
        }
      }
      {
        integer NR{m_nr};
        integer NC{m_qx};
        integer NN{NR*NC};
        for( integer i{0}; i < NN; ++i ) {
          integer r{ i % NR };
          integer c{ integer(i / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_Cqmat[i]),
            "BorderedCR::check_matrix, found {} at position ({},{}) of {}x{} matrix Cq\n",
            m_Cqmat[i], r, c, NR, NC
         );
        }
      }
      {
        integer NR{n};
        integer NC{n};
        integer SB{NR*NC};
        integer NN{SB*nblock};
        for( integer i{0}; i < NN; ++i ) {
          integer nb{ integer(i/SB) };
          integer r{ i % NR };
          integer c{ integer((i-nb*SB) / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_Dmat[i]),
            "BorderedCR::check_matrix, found {} at block {} position ({},{}) of {}x{} matrix D\n",
            m_Dmat[i], nb, r, c, NR, NC
         );
        }
      }
      {
        integer NR{n};
        integer NC{n};
        integer SB{NR*NC};
        integer NN{SB*nblock};
        for( integer i{0}; i < NN; ++i ) {
          integer nb{ integer(i/SB) };
          integer r{ i % NR };
          integer c{ integer((i-nb*SB) / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_Emat[i]),
            "BorderedCR::check_matrix, found {} at block {} position ({},{}) of {}x{} matrix E\n",
            m_Emat[i], nb, r, c, NR, NC
         );
        }
      }
      {
        integer NR{m_nr};
        integer NC{m_nx};
        integer NN{NR*NC};
        for( integer i{0}; i < NN; ++i ) {
          integer r{ i % NR };
          integer c{ integer(i / NR) };
          UTILS_ASSERT(
            Utils::is_finite(m_Fmat[i]),
            "BorderedCR::check_matrix, found {} at position ({},{}) of {}x{} matrix F\n",
            m_Fmat[i], r, c, NR, NC
         );
        }
      }
    } catch ( exception const & err ) {
      m_last_error = err.what();
      throw;
    } catch (...) {
      m_last_error = "BorderedCR::check_matrix, unknown error\n";
      throw;
    }
  }

  // ---------------------------------------------------------------------------

  template class BorderedCR<float>;
  template class BorderedCR<double>;

}
