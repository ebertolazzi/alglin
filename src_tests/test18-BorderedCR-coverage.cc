/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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
 |      Universita' degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

static Utils::Console msg( &std::cout, 4 );

using real_type = double;
using integer   = alglin::integer;
using BCR       = alglin::BorderedCR<real_type>;
using MatW      = alglin::MatrixWrapper<real_type>;
using Matrix    = alglin::Matrix<real_type>;

namespace {

  integer constexpr N_BLOCKS { 6 };
  integer constexpr N_SIZE   { 2 };
  integer constexpr QR_SIZE  { 1 };
  integer constexpr QX_SIZE  { 2 };
  integer constexpr NR_SIZE  { 3 };
  integer constexpr NX_SIZE  { 2 };

  real_type constexpr TOL{ 1e-10 };

  void
  check_close(
    real_type value,
    real_type expected,
    char const where[],
    integer i,
    integer j,
    real_type tol = TOL
  ) {
    UTILS_ASSERT(
      std::abs(value-expected) < tol,
      "{} mismatch at ({},{}): found {}, expected {}\n",
      where, i, j, value, expected
    );
  }

  void
  check_close(
    real_type value,
    real_type expected,
    char const where[],
    real_type tol = TOL
  ) {
    UTILS_ASSERT(
      std::abs(value-expected) < tol,
      "{} mismatch: found {}, expected {}\n",
      where, value, expected
    );
  }

  template <typename F>
  void
  expect_throw( F && fun, char const where[] ) {
    bool thrown{ false };
    try {
      fun();
    } catch ( ... ) {
      thrown = true;
    }
    UTILS_ASSERT( thrown, "{} expected an exception\n", where );
  }

  real_type
  tagged_value( integer tag, integer i, integer j ) {
    return real_type(tag) + real_type(10*j + i) / 10.0;
  }

  void
  fill_view( MatW & view, integer tag ) {
    for ( integer j{0}; j < view.ncols(); ++j )
      for ( integer i{0}; i < view.nrows(); ++i )
        view(i,j) = tagged_value( tag, i, j );
  }

  std::vector<real_type>
  dense_from_sparse( BCR const & bcr ) {
    integer const nnz{ bcr.sparse_nnz() };
    std::vector<integer>   I( static_cast<size_t>(nnz) );
    std::vector<integer>   J( static_cast<size_t>(nnz) );
    std::vector<real_type> V( static_cast<size_t>(nnz) );
    integer const nr{ bcr.nrows() };
    integer const nc{ bcr.ncols() };
    std::vector<real_type> A( static_cast<size_t>(nr * nc), 0 );

    bcr.sparse_pattern( I.data(), J.data(), 0 );
    bcr.sparse_values( V.data() );

    for ( integer k{0}; k < nnz; ++k ) {
      size_t const kk{ static_cast<size_t>(k) };
      A[ static_cast<size_t>( I[kk] + nr * J[kk] ) ] += V[kk];
    }

    return A;
  }

  std::vector<real_type>
  dense_mv(
    std::vector<real_type> const & A,
    integer                        nr,
    integer                        nc,
    real_type const                x[]
  ) {
    std::vector<real_type> y( static_cast<size_t>(nr), 0 );
    for ( integer j{0}; j < nc; ++j )
      for ( integer i{0}; i < nr; ++i )
        y[static_cast<size_t>(i)] += A[static_cast<size_t>(i + nr*j)] * x[j];
    return y;
  }

  void
  check_dense_equal(
    std::vector<real_type> const & A,
    std::vector<real_type> const & B,
    integer                        nr,
    integer                        nc,
    char const                     where[]
  ) {
    UTILS_ASSERT(
      A.size() == B.size(),
      "{} size mismatch: {} vs {}\n",
      where, A.size(), B.size()
    );
    for ( integer j{0}; j < nc; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close(
          A[static_cast<size_t>(i + nr*j)],
          B[static_cast<size_t>(i + nr*j)],
          where, i, j
        );
  }

  void
  fill_stable_system( BCR & bcr ) {
    integer const nblock{ bcr.number_of_blocks() };
    integer const n     { bcr.block_size()       };
    integer const qr    { bcr.dim_qr()           };
    integer const qx    { bcr.dim_qx()           };
    integer const nr    { bcr.dim_nr()           };
    integer const nx    { bcr.dim_nx()           };
    integer const m     { n + qr                 };

    bcr.fill_zero();

    for ( integer k{0}; k < nblock; ++k ) {
      bcr.D(k,0,0) = 5.0 + 0.2*k;
      bcr.D(k,1,0) = 0.15;
      bcr.D(k,0,1) = -0.08;
      bcr.D(k,1,1) = 4.7 + 0.2*k;

      bcr.E(k,0,0) = -0.35;
      bcr.E(k,1,0) = 0.04;
      bcr.E(k,0,1) = 0.06;
      bcr.E(k,1,1) = -0.28;

      for ( integer j{0}; j < nx; ++j )
        for ( integer i{0}; i < n; ++i )
          bcr.B(k,i,j) = 0.04 * (1 + k + i + 2*j);

      for ( integer j{0}; j < n; ++j )
        for ( integer i{0}; i < nr; ++i )
          bcr.C(k,i,j) = 0.03 * (1 + 2*k + i + j);
    }

    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i )
        bcr.C(nblock,i,j) = 0.05 * (1 + i + j);

    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        bcr.Cq(i,j) = 0.07 * (1 + i + 2*j);

    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        bcr.F(i,j) = 0.02 * (1 + i + j);

    bcr.F(0,0) += 3.2;
    bcr.F(1,1) += 2.8;

    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,j) = 0.02 * (1 + i + j);

    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,n+j) = 0.03 * (1 + i + j);

    bcr.H(0,n+0) += 5.1;
    bcr.H(1,n+1) += 4.9;
    bcr.H(2,n+0) += 0.2;
    bcr.H(2,n+1) -= 0.25;

    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,2*n+j) = 0.04 * (1 + i + j);

    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,2*n+qx+j) = 0.05 * (1 + i + j);
  }

  void
  build_square_solver_matrix( BCR & bcr ) {
    bcr.allocate( N_BLOCKS, N_SIZE, QR_SIZE, QX_SIZE, NR_SIZE, NX_SIZE );
    fill_stable_system( bcr );
  }

  void
  build_singular_fallback_matrix( BCR & bcr ) {
    bcr.allocate( 1, 1, 0, 0, 1, 1 );
    bcr.fill_zero();

    bcr.D(0,0,0) = 1;
    bcr.E(0,0,0) = 0;
    bcr.B(0,0,0) = 0;

    bcr.H(0,0) = 0;
    bcr.H(0,1) = 1;
    bcr.H(0,2) = 0;

    bcr.C(0,0,0) = 1;
    bcr.C(1,0,0) = 1;
    bcr.F(0,0)   = 0;
  }

  void
  verify_view( MatW const & view, integer rows, integer cols, char const where[] ) {
    UTILS_ASSERT(
      view.nrows() == rows && view.ncols() == cols,
      "{} wrong size {} x {}, expected {} x {}\n",
      where, view.nrows(), view.ncols(), rows, cols
    );
  }

  void
  test_api_and_selection() {
    Utils::ThreadPool1 TP(4);
    BCR                bcr( &TP, 1 );
    BCR                bcr_fake( nullptr );

    bcr.allocate( N_BLOCKS, N_SIZE, QR_SIZE, QX_SIZE, NR_SIZE, NX_SIZE );
    bcr_fake.allocate( 2, 2, 0, 0, 0, 0 );

    UTILS_ASSERT0( bcr.number_of_blocks() == N_BLOCKS, "wrong number_of_blocks\n" );
    UTILS_ASSERT0( bcr.block_size()       == N_SIZE,   "wrong block_size\n" );
    UTILS_ASSERT0( bcr.dim_qr()           == QR_SIZE,  "wrong qr size\n" );
    UTILS_ASSERT0( bcr.dim_qx()           == QX_SIZE,  "wrong qx size\n" );
    UTILS_ASSERT0( bcr.dim_nr()           == NR_SIZE,  "wrong nr size\n" );
    UTILS_ASSERT0( bcr.dim_nx()           == NX_SIZE,  "wrong nx size\n" );
    UTILS_ASSERT0( bcr.Nr() == 2 * N_SIZE + QR_SIZE + NR_SIZE, "reduced Nr mismatch\n" );
    UTILS_ASSERT0( bcr.Nc() == 2 * N_SIZE + QX_SIZE + NX_SIZE, "reduced Nc mismatch\n" );
    UTILS_ASSERT0( bcr.nrows() == N_SIZE * (N_BLOCKS + 1) + QR_SIZE + NR_SIZE, "global nrows mismatch\n" );
    UTILS_ASSERT0( bcr.ncols() == N_SIZE * (N_BLOCKS + 1) + QX_SIZE + NX_SIZE, "global ncols mismatch\n" );

    bcr.set_factorize_use_thread( false );
    bcr.set_solve_use_thread( false );
    UTILS_ASSERT0( !bcr.factorize_use_thread(), "factorize thread flag mismatch\n" );
    UTILS_ASSERT0( !bcr.solve_use_thread(),     "solve thread flag mismatch\n" );
    bcr.set_factorize_use_thread( true );
    bcr.set_solve_use_thread( true );
    UTILS_ASSERT0( bcr.factorize_use_thread(), "factorize thread flag mismatch\n" );
    UTILS_ASSERT0( bcr.solve_use_thread(),     "solve thread flag mismatch\n" );

    bcr.select( "LU" );
    UTILS_ASSERT0( bcr.info_algo_block() == "CyclicReduction+LU", "select(LU) failed\n" );
    bcr.select( "QR" );
    UTILS_ASSERT0( bcr.info_algo_block() == "CyclicReduction+QR", "select(QR) failed\n" );
    bcr.select_QRP();
    UTILS_ASSERT0( bcr.info_algo_block() == "CyclicReduction+QRP", "select_QRP failed\n" );

    bcr.select_last( "LU" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock LU", "select_last(LU) failed\n" );
    bcr.select_last( "LUPQ" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock LUPQ", "select_last(LUPQ) failed\n" );
    bcr.select_last( "QR" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock QR", "select_last(QR) failed\n" );
    bcr.select_last( "QRP" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock QRP", "select_last(QRP) failed\n" );
    bcr.select_last( "SVD" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock SVD", "select_last(SVD) failed\n" );
    bcr.select_last( "LSS" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock LSS", "select_last(LSS) failed\n" );
    bcr.select_last( "LSY" );
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock LSY", "select_last(LSY) failed\n" );
    bcr.select_last_PINV();
    UTILS_ASSERT0( bcr.info_algo_last_block() == "LastBlock PINV", "select_last_PINV failed\n" );

    bcr.select_last2( "NONE" );
    UTILS_ASSERT0( bcr.info_algo_last_block2() == "LastBlock NONE", "select_last2(NONE) failed\n" );
    bcr.select_last2( "SVD" );
    UTILS_ASSERT0( bcr.info_algo_last_block2() == "LastBlock SVD", "select_last2(SVD) failed\n" );
    bcr.select_last2( "LSS" );
    UTILS_ASSERT0( bcr.info_algo_last_block2() == "LastBlock LSS", "select_last2(LSS) failed\n" );
    bcr.select_last2( "LSY" );
    UTILS_ASSERT0( bcr.info_algo_last_block2() == "LastBlock LSY", "select_last2(LSY) failed\n" );
    bcr.select_last2_PINV();
    UTILS_ASSERT0( bcr.info_algo_last_block2() == "LastBlock PINV", "select_last2_PINV failed\n" );

    UTILS_ASSERT0(
      BCR::choice_to_string( BCR::BORDERED_Choice::LU ) == "CyclicReduction+LU",
      "choice_to_string(BORDERED_Choice::LU) failed\n"
    );
    UTILS_ASSERT0(
      BCR::choice_to_string( BCR::BORDERED_LAST_Choice::QRP ) == "LastBlock QRP",
      "choice_to_string(BORDERED_LAST_Choice::QRP) failed\n"
    );
    UTILS_ASSERT0(
      BCR::choice_to_string( BCR::BORDERED_LAST_Choice2::NONE ) == "LastBlock NONE",
      "choice_to_string(BORDERED_LAST_Choice2::NONE) failed\n"
    );

    bcr.set_num_parallel_block( 4 );
    std::string const info{ bcr.info("  ") };
    UTILS_ASSERT0( info.find("nblk:3") != std::string::npos, "set_num_parallel_block not reflected in info()\n" );
    UTILS_ASSERT0( info.find("CyclicReduction+QRP") != std::string::npos, "info() missing internal solver\n" );
    UTILS_ASSERT0( info.find("LastBlock PINV") != std::string::npos, "info() missing last solver\n" );

    std::ostringstream sout;
    bcr.info( sout );
    UTILS_ASSERT0( sout.str().find("rows") != std::string::npos, "info(stream) produced empty output\n" );

    expect_throw( [&]() { bcr.select("BAD"); },        "select invalid option" );
    expect_throw( [&]() { bcr.select_last("BAD"); },   "select_last invalid option" );
    expect_throw( [&]() { bcr.select_last2("BAD"); },  "select_last2 invalid option" );
  }

  void
  test_zero_and_load_methods() {
    Utils::ThreadPool1 TP(4);
    BCR                bcr( &TP );

    integer const nblock{ 3 };
    integer const n     { 2 };
    integer const qr    { 1 };
    integer const qx    { 2 };
    integer const nr    { 3 };
    integer const nx    { 2 };
    integer const m     { n + qr };

    bcr.allocate( nblock, n, qr, qx, nr, nx );
    bcr.fill_zero();

    bcr.D(0,0,0) = 1;
    bcr.zero_D();
    check_close( bcr.D(0,0,0), 0.0, "zero_D" );

    bcr.E(0,0,0) = 2;
    bcr.zero_E();
    check_close( bcr.E(0,0,0), 0.0, "zero_E" );

    bcr.B(0,0,0) = 3;
    bcr.zero_B();
    check_close( bcr.B(0,0,0), 0.0, "zero_B" );

    bcr.F(0,0) = 4;
    bcr.zero_F();
    check_close( bcr.F(0,0), 0.0, "zero_F" );

    bcr.H(0,0) = 5;
    bcr.zero_H();
    check_close( bcr.H(0,0), 0.0, "zero_H" );

    bcr.C(0,0,0) = 6;
    bcr.zero_C();
    check_close( bcr.C(0,0,0), 0.0, "zero_C" );

    bcr.Cq(0,0) = 7;
    bcr.zero_Cq();
    check_close( bcr.Cq(0,0), 0.0, "zero_Cq" );

    bcr.D(0,0,0)  = 1;
    bcr.E(0,0,0)  = 1;
    bcr.B(0,0,0)  = 1;
    bcr.F(0,0)    = 1;
    bcr.H(0,0)    = 1;
    bcr.C(0,0,0)  = 1;
    bcr.Cq(0,0)   = 1;
    bcr.fill_zero();
    check_close( bcr.D(0,0,0), 0.0, "fill_zero/D" );
    check_close( bcr.E(0,0,0), 0.0, "fill_zero/E" );
    check_close( bcr.B(0,0,0), 0.0, "fill_zero/B" );
    check_close( bcr.F(0,0),   0.0, "fill_zero/F" );
    check_close( bcr.H(0,0),   0.0, "fill_zero/H" );
    check_close( bcr.C(0,0,0), 0.0, "fill_zero/C" );
    check_close( bcr.Cq(0,0),  0.0, "fill_zero/Cq" );

    Matrix B_storage( n+2, nx );
    MatW   B_view;
    B_storage.fill( -1 );
    B_storage.view_block( 1, 0, n, nx, B_view );
    fill_view( B_view, 10 );
    bcr.load_B( 1, B_view );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.B(1,i,j), B_view(i,j), "load_B(MatW)", i, j );

    Matrix B_add_storage( n+1, nx );
    MatW   B_add_view;
    B_add_storage.fill( -2 );
    B_add_storage.view_block( 0, 0, n, nx, B_add_view );
    fill_view( B_add_view, 20 );
    bcr.add_to_B( 1, B_add_view.data(), B_add_view.ldim() );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.B(1,i,j), B_view(i,j) + B_add_view(i,j), "add_to_B(raw)", i, j );

    Matrix B_add2_storage( n+1, nx );
    MatW   B_add2_view;
    B_add2_storage.view_block( 0, 0, n, nx, B_add2_view );
    fill_view( B_add2_view, 30 );
    bcr.add_to_B( 1, B_add2_view );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close(
          bcr.B(1,i,j),
          B_view(i,j) + B_add_view(i,j) + B_add2_view(i,j),
          "add_to_B(MatW)", i, j
        );

    Matrix C_storage( nr+1, n );
    MatW   C_view;
    C_storage.view_block( 0, 0, nr, n, C_view );
    fill_view( C_view, 40 );
    bcr.load_C( 2, C_view.data(), C_view.ldim() );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.C(2,i,j), C_view(i,j), "load_C(raw)", i, j );

    Matrix C_add_storage( nr+2, n );
    MatW   C_add_view;
    C_add_storage.view_block( 1, 0, nr, n, C_add_view );
    fill_view( C_add_view, 50 );
    bcr.add_to_C( 2, C_add_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.C(2,i,j), C_view(i,j) + C_add_view(i,j), "add_to_C(MatW)", i, j );

    Matrix C_add2_storage( nr+1, n );
    MatW   C_add2_view;
    C_add2_storage.view_block( 0, 0, nr, n, C_add2_view );
    fill_view( C_add2_view, 60 );
    bcr.add_to_C( 2, C_add2_view.data(), C_add2_view.ldim() );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close(
          bcr.C(2,i,j),
          C_view(i,j) + C_add_view(i,j) + C_add2_view(i,j),
          "add_to_C(raw)", i, j
        );

    Matrix D_storage( n+1, n );
    MatW   D_view;
    D_storage.view_block( 0, 0, n, n, D_view );
    fill_view( D_view, 70 );
    bcr.load_D( 0, D_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.D(0,i,j), D_view(i,j), "load_D(MatW)", i, j );

    Matrix D_raw_storage( n+2, n );
    MatW   D_raw_view;
    D_raw_storage.view_block( 1, 0, n, n, D_raw_view );
    fill_view( D_raw_view, 80 );
    bcr.load_D( 1, D_raw_view.data(), D_raw_view.ldim() );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.D(1,i,j), D_raw_view(i,j), "load_D(raw)", i, j );

    Matrix E_storage( n+2, n );
    MatW   E_view;
    E_storage.view_block( 1, 0, n, n, E_view );
    fill_view( E_view, 90 );
    bcr.load_E( 0, E_view.data(), E_view.ldim() );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.E(0,i,j), E_view(i,j), "load_E(raw)", i, j );

    Matrix E_raw_storage( n+1, n );
    MatW   E_raw_view;
    E_raw_storage.view_block( 0, 0, n, n, E_raw_view );
    fill_view( E_raw_view, 100 );
    bcr.load_E( 1, E_raw_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.E(1,i,j), E_raw_view(i,j), "load_E(MatW)", i, j );

    Matrix DE_storage( n+1, 2*n );
    MatW   DE_view;
    DE_storage.view_block( 0, 0, n, 2*n, DE_view );
    fill_view( DE_view, 110 );
    bcr.load_DE( 2, DE_view.data(), DE_view.ldim() );
    for ( integer j{0}; j < n; ++j ) {
      for ( integer i{0}; i < n; ++i ) {
        check_close( bcr.D(2,i,j), DE_view(i,j),    "load_DE/D", i, j );
        check_close( bcr.E(2,i,j), DE_view(i,n+j),  "load_DE/E", i, j );
      }
    }

    Matrix DEB_storage( n+2, 2*n+nx );
    MatW   DEB_view;
    DEB_storage.view_block( 1, 0, n, 2*n+nx, DEB_view );
    fill_view( DEB_view, 120 );
    bcr.load_DEB( 1, DEB_view.data(), DEB_view.ldim() );
    for ( integer j{0}; j < n; ++j ) {
      for ( integer i{0}; i < n; ++i ) {
        check_close( bcr.D(1,i,j), DEB_view(i,j),   "load_DEB/D", i, j );
        check_close( bcr.E(1,i,j), DEB_view(i,n+j), "load_DEB/E", i, j );
      }
    }
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( bcr.B(1,i,j), DEB_view(i,2*n+j), "load_DEB/B", i, j );

    Matrix F_storage( nr+2, nx );
    MatW   F_view;
    F_storage.view_block( 1, 0, nr, nx, F_view );
    fill_view( F_view, 130 );
    bcr.load_F( F_view.data(), F_view.ldim() );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j), F_view(i,j), "load_F(raw)", i, j );

    Matrix F_add_storage( nr+1, nx );
    MatW   F_add_view;
    F_add_storage.view_block( 0, 0, nr, nx, F_add_view );
    fill_view( F_add_view, 140 );
    bcr.add_to_F( F_add_view );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j), F_view(i,j) + F_add_view(i,j), "add_to_F(MatW)", i, j );

    Matrix F_add2_storage( nr+1, nx );
    MatW   F_add2_view;
    F_add2_storage.view_block( 0, 0, nr, nx, F_add2_view );
    fill_view( F_add2_view, 150 );
    bcr.add_to_F( F_add2_view.data(), F_add2_view.ldim() );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close(
          bcr.F(i,j),
          F_view(i,j) + F_add_view(i,j) + F_add2_view(i,j),
          "add_to_F(raw)", i, j
        );

    Matrix Cq_storage( nr+1, qx );
    MatW   Cq_view;
    Cq_storage.view_block( 0, 0, nr, qx, Cq_view );
    fill_view( Cq_view, 160 );
    bcr.load_Cq( Cq_view );
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.Cq(i,j), Cq_view(i,j), "load_Cq(MatW)", i, j );

    Matrix Cq_raw_storage( nr+2, qx );
    MatW   Cq_raw_view;
    Cq_raw_storage.view_block( 1, 0, nr, qx, Cq_raw_view );
    fill_view( Cq_raw_view, 170 );
    bcr.load_Cq( Cq_raw_view.data(), Cq_raw_view.ldim() );
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.Cq(i,j), Cq_raw_view(i,j), "load_Cq(raw)", i, j );

    Matrix CqF_storage( nr+2, qx+nx );
    MatW   CqF_view;
    CqF_storage.view_block( 1, 0, nr, qx+nx, CqF_view );
    fill_view( CqF_view, 180 );
    bcr.load_CqF( CqF_view.data(), CqF_view.ldim() );
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.Cq(i,j), CqF_view(i,j), "load_CqF/Cq", i, j );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j), CqF_view(i,qx+j), "load_CqF/F", i, j );

    Matrix H0_storage( m+1, n );
    Matrix HN_storage( m+2, n );
    Matrix Hq_storage( m+1, qx );
    Matrix Hp_storage( m+1, nx );
    MatW   H0_view, HN_view, Hq_view, Hp_view;
    H0_storage.view_block( 0, 0, m, n,  H0_view );
    HN_storage.view_block( 1, 0, m, n,  HN_view );
    Hq_storage.view_block( 0, 0, m, qx, Hq_view );
    Hp_storage.view_block( 0, 0, m, nx, Hp_view );
    fill_view( H0_view, 190 );
    fill_view( HN_view, 200 );
    fill_view( Hq_view, 210 );
    fill_view( Hp_view, 220 );
    bcr.load_bottom(
      H0_view.data(), H0_view.ldim(),
      HN_view.data(), HN_view.ldim(),
      Hq_view.data(), Hq_view.ldim(),
      Hp_view.data(), Hp_view.ldim()
    );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < m; ++i ) {
        check_close( bcr.H(i,j),          H0_view(i,j), "load_bottom(raw)/H0", i, j );
        check_close( bcr.H(i,n+j),        HN_view(i,j), "load_bottom(raw)/HN", i, j );
      }
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < m; ++i )
        check_close( bcr.H(i,2*n+j),      Hq_view(i,j), "load_bottom(raw)/Hq", i, j );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < m; ++i )
        check_close( bcr.H(i,2*n+qx+j),   Hp_view(i,j), "load_bottom(raw)/Hp", i, j );

    fill_view( H0_view, 230 );
    fill_view( HN_view, 240 );
    fill_view( Hq_view, 250 );
    fill_view( Hp_view, 260 );
    bcr.load_bottom( H0_view, HN_view, Hq_view, Hp_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < m; ++i ) {
        check_close( bcr.H(i,j),        H0_view(i,j), "load_bottom(MatW4)/H0", i, j );
        check_close( bcr.H(i,n+j),      HN_view(i,j), "load_bottom(MatW4)/HN", i, j );
      }

    Matrix H_storage( m+2, 2*n+qx+nx );
    MatW   H_view;
    H_storage.view_block( 1, 0, m, 2*n+qx+nx, H_view );
    fill_view( H_view, 270 );
    bcr.load_bottom( H_view.data(), H_view.ldim() );
    for ( integer j{0}; j < H_view.ncols(); ++j )
      for ( integer i{0}; i < H_view.nrows(); ++i )
        check_close( bcr.H(i,j), H_view(i,j), "load_bottom(rawH)", i, j );

    fill_view( H_view, 280 );
    bcr.load_bottom( H_view );
    for ( integer j{0}; j < H_view.ncols(); ++j )
      for ( integer i{0}; i < H_view.nrows(); ++i )
        check_close( bcr.H(i,j), H_view(i,j), "load_bottom(MatWH)", i, j );

    Matrix C0_storage( nr+2, n );
    Matrix CN_storage( nr+2, n );
    Matrix CBq_storage( nr+1, qx );
    Matrix CF_storage( nr+1, nx );
    MatW   C0_view, CN_view, Cqb_view, Fb_view;
    C0_storage.view_block( 1, 0, nr, n,  C0_view );
    CN_storage.view_block( 1, 0, nr, n,  CN_view );
    CBq_storage.view_block( 0, 0, nr, qx, Cqb_view );
    CF_storage.view_block( 0, 0, nr, nx, Fb_view );
    fill_view( C0_view, 290 );
    fill_view( CN_view, 300 );
    fill_view( Cqb_view, 310 );
    fill_view( Fb_view, 320 );
    bcr.load_bottom2(
      C0_view.data(), C0_view.ldim(),
      CN_view.data(), CN_view.ldim(),
      Cqb_view.data(), Cqb_view.ldim(),
      Fb_view.data(), Fb_view.ldim()
    );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(0,i,j),      C0_view(i,j), "load_bottom2(raw)/C0", i, j );
        check_close( bcr.C(nblock,i,j), CN_view(i,j), "load_bottom2(raw)/CN", i, j );
      }
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.Cq(i,j),       Cqb_view(i,j), "load_bottom2(raw)/Cq", i, j );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j),        Fb_view(i,j), "load_bottom2(raw)/F", i, j );

    fill_view( C0_view, 330 );
    fill_view( CN_view, 340 );
    fill_view( Cqb_view, 350 );
    fill_view( Fb_view, 360 );
    bcr.load_bottom2( C0_view, CN_view, Cqb_view, Fb_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(0,i,j),      C0_view(i,j), "load_bottom2(MatW4)/C0", i, j );
        check_close( bcr.C(nblock,i,j), CN_view(i,j), "load_bottom2(MatW4)/CN", i, j );
      }

    Matrix bottom2_storage( nr+2, 2*n+qx+nx );
    MatW   bottom2_view;
    bottom2_storage.view_block( 1, 0, nr, 2*n+qx+nx, bottom2_view );
    fill_view( bottom2_view, 370 );
    bcr.load_bottom2( bottom2_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(0,i,j),      bottom2_view(i,j),       "load_bottom2(MatWH)/C0", i, j );
        check_close( bcr.C(nblock,i,j), bottom2_view(i,n+j),     "load_bottom2(MatWH)/CN", i, j );
      }
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.Cq(i,j),       bottom2_view(i,2*n+j),   "load_bottom2(MatWH)/Cq", i, j );
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j),        bottom2_view(i,2*n+qx+j),"load_bottom2(MatWH)/F", i, j );

    bcr.zero_C();
    Matrix C2_storage( nr+1, 2*n );
    MatW   C2_view;
    C2_storage.view_block( 0, 0, nr, 2*n, C2_view );
    fill_view( C2_view, 380 );
    bcr.add_to_C2( 1, C2_view.data(), C2_view.ldim() );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(1,i,j), C2_view(i,j),   "add_to_C2(raw)/C1", i, j );
        check_close( bcr.C(2,i,j), C2_view(i,n+j), "add_to_C2(raw)/C2", i, j );
      }

    Matrix C2_add_storage( nr+2, 2*n );
    MatW   C2_add_view;
    C2_add_storage.view_block( 1, 0, nr, 2*n, C2_add_view );
    fill_view( C2_add_view, 390 );
    bcr.add_to_C2( 1, C2_add_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(1,i,j), C2_view(i,j)   + C2_add_view(i,j),   "add_to_C2(MatW)/C1", i, j );
        check_close( bcr.C(2,i,j), C2_view(i,n+j) + C2_add_view(i,n+j), "add_to_C2(MatW)/C2", i, j );
      }

    bcr.zero_C();
    bcr.zero_F();
    Matrix C2F_storage( nr+1, 2*n+nx );
    MatW   C2F_view;
    C2F_storage.view_block( 0, 0, nr, 2*n+nx, C2F_view );
    fill_view( C2F_view, 400 );
    bcr.add_to_C2F( 1, C2F_view.data(), C2F_view.ldim() );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(1,i,j), C2F_view(i,j),   "add_to_C2F(raw)/C1", i, j );
        check_close( bcr.C(2,i,j), C2F_view(i,n+j), "add_to_C2F(raw)/C2", i, j );
      }
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j), C2F_view(i,2*n+j), "add_to_C2F(raw)/F", i, j );

    Matrix C2F_add_storage( nr+2, 2*n+nx );
    MatW   C2F_add_view;
    C2F_add_storage.view_block( 1, 0, nr, 2*n+nx, C2F_add_view );
    fill_view( C2F_add_view, 410 );
    bcr.add_to_C2F( 1, C2F_add_view );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i ) {
        check_close( bcr.C(1,i,j), C2F_view(i,j)   + C2F_add_view(i,j),   "add_to_C2F(MatW)/C1", i, j );
        check_close( bcr.C(2,i,j), C2F_view(i,n+j) + C2F_add_view(i,n+j), "add_to_C2F(MatW)/C2", i, j );
      }
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( bcr.F(i,j), C2F_view(i,2*n+j) + C2F_add_view(i,2*n+j), "add_to_C2F(MatW)/F", i, j );

    MatW B_wrap, C_wrap, D_wrap, E_wrap, F_wrap, Cq_wrap, H_wrap;
    bcr.B( 1, B_wrap );
    bcr.C( 1, C_wrap );
    bcr.D( 1, D_wrap );
    bcr.E( 1, E_wrap );
    bcr.F(    F_wrap );
    bcr.Cq(   Cq_wrap );
    bcr.H(    H_wrap );

    verify_view( B_wrap, n,  nx, "B(MatW)" );
    verify_view( C_wrap, nr, n,  "C(MatW)" );
    verify_view( D_wrap, n,  n,  "D(MatW)" );
    verify_view( E_wrap, n,  n,  "E(MatW)" );
    verify_view( F_wrap, nr, nx, "F(MatW)" );
    verify_view( Cq_wrap, nr, qx, "Cq(MatW)" );
    verify_view( H_wrap, m, 2*n+qx+nx, "H(MatW)" );

    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < n; ++i )
        check_close( B_wrap(i,j), bcr.B(1,i,j), "B(MatW) values", i, j );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( C_wrap(i,j), bcr.C(1,i,j), "C(MatW) values", i, j );
    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < n; ++i ) {
        check_close( D_wrap(i,j), bcr.D(1,i,j), "D(MatW) values", i, j );
        check_close( E_wrap(i,j), bcr.E(1,i,j), "E(MatW) values", i, j );
      }
    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( F_wrap(i,j), bcr.F(i,j), "F(MatW) values", i, j );
    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        check_close( Cq_wrap(i,j), bcr.Cq(i,j), "Cq(MatW) values", i, j );
    for ( integer j{0}; j < H_wrap.ncols(); ++j )
      for ( integer i{0}; i < H_wrap.nrows(); ++i )
        check_close( H_wrap(i,j), bcr.H(i,j), "H(MatW) values", i, j );
  }

  void
  test_sparse_mv_dup_and_print() {
    Utils::ThreadPool1 TP(4);
    BCR                bcr( &TP );
    BCR                dup( &TP );
    BCR                from_sparse( &TP );

    build_square_solver_matrix( bcr );
    bcr.check_matrix();

    dup.dup( bcr );
    check_dense_equal(
      dense_from_sparse( bcr ),
      dense_from_sparse( dup ),
      bcr.nrows(),
      bcr.ncols(),
      "dup sparse equality"
    );

    integer const nnz{ bcr.sparse_nnz() };
    std::vector<integer>   I( static_cast<size_t>(nnz) );
    std::vector<integer>   J( static_cast<size_t>(nnz) );
    std::vector<real_type> V( static_cast<size_t>(nnz) );
    bcr.sparse_pattern( I.data(), J.data(), 0 );
    bcr.sparse_values( V.data() );

    from_sparse.allocate(
      bcr.number_of_blocks(),
      bcr.block_size(),
      bcr.dim_qr(),
      bcr.dim_qx(),
      bcr.dim_nr(),
      bcr.dim_nx()
    );
    from_sparse.sparse_load( V.data(), I.data(), 0, J.data(), 0, nnz );

    check_dense_equal(
      dense_from_sparse( bcr ),
      dense_from_sparse( from_sparse ),
      bcr.nrows(),
      bcr.ncols(),
      "sparse_load roundtrip"
    );

    std::vector<real_type> const A{ dense_from_sparse( bcr ) };
    integer const nr{ bcr.nrows() };
    integer const nc{ bcr.ncols() };

    std::vector<real_type>       x( static_cast<size_t>(nc) );
    std::vector<real_type>       y( static_cast<size_t>(nr) );
    std::vector<real_type>       y_dense;

    for ( integer i{0}; i < nc; ++i ) x[static_cast<size_t>(i)] = 1.0 + 0.25*i;

    bcr.Mv( x.data(), y.data() );
    y_dense = dense_mv( A, nr, nc, x.data() );
    for ( integer i{0}; i < nr; ++i )
      check_close( y[static_cast<size_t>(i)], y_dense[static_cast<size_t>(i)], "Mv", i, 0 );

    std::vector<real_type> acc( static_cast<size_t>(nr), 3.0 );
    std::vector<real_type> acc_dense( static_cast<size_t>(nr), 3.0 );
    bcr.add_Mv( 0.5, x.data(), acc.data() );
    y_dense = dense_mv( A, nr, nc, x.data() );
    for ( integer i{0}; i < nr; ++i ) acc_dense[static_cast<size_t>(i)] += 0.5 * y_dense[static_cast<size_t>(i)];
    for ( integer i{0}; i < nr; ++i )
      check_close( acc[static_cast<size_t>(i)], acc_dense[static_cast<size_t>(i)], "add_Mv(alpha)", i, 0 );

    std::vector<real_type> acc2( static_cast<size_t>(nr), -2.0 );
    std::vector<real_type> acc2_dense( static_cast<size_t>(nr), -2.0 );
    bcr.add_Mv( x.data(), acc2.data() );
    for ( integer i{0}; i < nr; ++i ) acc2_dense[static_cast<size_t>(i)] += y_dense[static_cast<size_t>(i)];
    for ( integer i{0}; i < nr; ++i )
      check_close( acc2[static_cast<size_t>(i)], acc2_dense[static_cast<size_t>(i)], "add_Mv()", i, 0 );

    std::ostringstream sout;
    bcr.print_matlab_script( sout );
    std::string const matlab{ sout.str() };
    UTILS_ASSERT0( matlab.find("MAT = sparse") != std::string::npos, "print_matlab_script missing sparse() call\n" );
    UTILS_ASSERT0( matlab.find("nr  = ") != std::string::npos, "print_matlab_script missing nr\n" );
    UTILS_ASSERT0( matlab.find("nc  = ") != std::string::npos, "print_matlab_script missing nc\n" );
  }

  void
  test_check_matrix_and_tsolve() {
    Utils::ThreadPool1 TP(4);
    BCR                ok( &TP );
    BCR                bad( &TP );

    build_square_solver_matrix( ok );
    ok.set_debug( true );
    ok.check_matrix();
    ok.set_debug( false );

    build_square_solver_matrix( bad );
    bad.F( 0, 0 ) = std::numeric_limits<real_type>::quiet_NaN();
    expect_throw( [&]() { bad.check_matrix(); }, "check_matrix on NaN" );
    UTILS_ASSERT0( bad.last_error().find("check_matrix") != std::string::npos, "last_error not updated after check_matrix failure\n" );

    std::vector<real_type> x( static_cast<size_t>(ok.ncols()), 0 );
    expect_throw( [&]() { ok.t_solve( x.data() ); }, "t_solve(single)" );
    expect_throw( [&]() { ok.t_solve( 2, x.data(), ok.ncols() ); }, "t_solve(multi)" );

    Matrix X( std::max(ok.nrows(), ok.ncols()), 2 );
    X.fill( 0 );
    expect_throw( [&]() { ok.t_solve( X ); }, "t_solve(MatrixWrapper)" );
  }

  void
  run_single_and_multi_solve_case(
    char const internal_choice[],
    char const last_choice[],
    bool       factorize_threads,
    bool       solve_threads
  ) {
    Utils::ThreadPool1 TP(4);
    BCR                rhs_builder( &TP );
    BCR                solver( &TP );

    build_square_solver_matrix( rhs_builder );
    build_square_solver_matrix( solver );

    solver.select( internal_choice );
    solver.select_last( last_choice );
    solver.set_num_parallel_block( 4 );
    solver.set_factorize_use_thread( factorize_threads );
    solver.set_solve_use_thread( solve_threads );
    solver.set_debug( true );

    integer const Nr{ solver.nrows() };
    integer const Nc{ solver.ncols() };
    integer const N { std::max(Nr,Nc) };
    integer const nrhs{ 3 };

    std::vector<real_type> xref( static_cast<size_t>(Nc) );
    std::vector<real_type> rhs( static_cast<size_t>(Nr) );
    std::vector<real_type> x( static_cast<size_t>(N) );
    std::vector<real_type> resid( static_cast<size_t>(Nr) );
    std::vector<real_type> Xref( static_cast<size_t>(Nc * nrhs) );
    std::vector<real_type> RHS( static_cast<size_t>(Nr * nrhs) );
    std::vector<real_type> X( static_cast<size_t>(N * nrhs) );

    for ( integer i{0}; i < Nc; ++i )
      xref[static_cast<size_t>(i)] = 0.25 + 0.5 * i;

    rhs_builder.Mv( xref.data(), rhs.data() );

    bool const ok_factorize{ solver.factorize() };
    UTILS_ASSERT(
      ok_factorize,
      "factorize failed for internal={}, last={}, fthread={}, sthread={}: {}\n",
      internal_choice, last_choice, factorize_threads, solve_threads, solver.last_error()
    );

    x = rhs;
    bool const ok_solve{ solver.solve( x.data() ) };
    UTILS_ASSERT(
      ok_solve,
      "solve failed for internal={}, last={}, fthread={}, sthread={}: {}\n",
      internal_choice, last_choice, factorize_threads, solve_threads, solver.last_error()
    );
    for ( integer i{0}; i < Nc; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "solve(single)", i, 0, 5e-9 );

    rhs_builder.Mv( x.data(), resid.data() );
    for ( integer i{0}; i < Nr; ++i )
      check_close( resid[static_cast<size_t>(i)], rhs[static_cast<size_t>(i)], "solve(single) residual", i, 0, 5e-9 );

    for ( integer c{0}; c < nrhs; ++c ) {
      for ( integer i{0}; i < Nc; ++i )
        Xref[static_cast<size_t>( i + c * Nc )] = (1 + c) * xref[static_cast<size_t>(i)] - 0.1 * c;
      rhs_builder.Mv( Xref.data() + c * Nc, RHS.data() + c * Nr );
    }

    for ( integer c{0}; c < nrhs; ++c )
      for ( integer i{0}; i < Nr; ++i )
        X[static_cast<size_t>( i + c * N )] = RHS[static_cast<size_t>( i + c * Nr )];
    bool const ok_solve_n{ solver.solve( nrhs, X.data(), N ) };
    UTILS_ASSERT(
      ok_solve_n,
      "solve(nrhs) failed for internal={}, last={}, fthread={}, sthread={}: {}\n",
      internal_choice, last_choice, factorize_threads, solve_threads, solver.last_error()
    );

    for ( integer c{0}; c < nrhs; ++c ) {
      for ( integer i{0}; i < Nc; ++i )
        check_close(
          X[static_cast<size_t>( i + c * N )],
          Xref[static_cast<size_t>( i + c * Nc )],
          "solve(nrhs)", i, c, 5e-9
        );
      rhs_builder.Mv( X.data() + c * N, resid.data() );
      for ( integer i{0}; i < Nr; ++i )
        check_close(
          resid[static_cast<size_t>(i)],
          RHS[static_cast<size_t>( i + c * Nr )],
          "solve(nrhs) residual", i, c, 5e-9
        );
    }

    Matrix Xmat( N, nrhs );
    Xmat.fill( 0 );
    for ( integer c{0}; c < nrhs; ++c )
      for ( integer i{0}; i < Nr; ++i )
        Xmat(i,c) = RHS[static_cast<size_t>( i + c * Nr )];

    bool const ok_solve_mat{ solver.solve( Xmat ) };
    UTILS_ASSERT(
      ok_solve_mat,
      "solve(MatrixWrapper) failed for internal={}, last={}, fthread={}, sthread={}: {}\n",
      internal_choice, last_choice, factorize_threads, solve_threads, solver.last_error()
    );

    for ( integer c{0}; c < nrhs; ++c ) {
      for ( integer i{0}; i < Nc; ++i )
        check_close(
          Xmat(i,c),
          Xref[static_cast<size_t>( i + c * Nc )],
          "solve(MatrixWrapper)", i, c, 5e-9
        );
      rhs_builder.Mv( Xmat.data() + c * Xmat.ldim(), resid.data() );
      for ( integer i{0}; i < Nr; ++i )
        check_close(
          resid[static_cast<size_t>(i)],
          RHS[static_cast<size_t>( i + c * Nr )],
          "solve(MatrixWrapper) residual", i, c, 5e-9
        );
    }
  }

  void
  test_solver_combinations() {
    std::array<char const *,3> const internals{
      "LU", "QR", "QRP"
    };
    std::array<char const *,8> const lasts{
      "LU", "LUPQ", "QR", "QRP", "SVD", "LSS", "LSY", "PINV"
    };

    for ( auto internal_choice : internals ) {
      for ( auto last_choice : lasts ) {
        run_single_and_multi_solve_case( internal_choice, last_choice, false, false );
        run_single_and_multi_solve_case( internal_choice, last_choice, true,  false );
        run_single_and_multi_solve_case( internal_choice, last_choice, false, true  );
        run_single_and_multi_solve_case( internal_choice, last_choice, true,  true  );
      }
    }
  }

  void
  test_last_solver_fallback() {
    Utils::ThreadPool1 TP(4);
    BCR                base( &TP );
    BCR                no_fallback( &TP );
    BCR                with_fallback( &TP );

    build_singular_fallback_matrix( base );
    build_singular_fallback_matrix( no_fallback );
    build_singular_fallback_matrix( with_fallback );

    no_fallback.select_LU();
    no_fallback.select_last_LU();
    no_fallback.select_last2_NONE();
    UTILS_ASSERT0( !no_fallback.factorize(), "singular LU case should fail without fallback\n" );
    UTILS_ASSERT0( no_fallback.last_error().find("load_and_factorize_last") != std::string::npos, "unexpected last_error without fallback\n" );

    with_fallback.select_LU();
    with_fallback.select_last_LU();
    with_fallback.select_last2_PINV();
    UTILS_ASSERT0( with_fallback.factorize(), "PINV fallback should recover singular last block\n" );

    std::vector<real_type> xref{ 2.0, -1.0, 0.0 };
    std::vector<real_type> rhs( 3 );
    std::vector<real_type> x( 3 );
    std::vector<real_type> resid( 3 );

    base.Mv( xref.data(), rhs.data() );
    x = rhs;

    UTILS_ASSERT0( with_fallback.solve( x.data() ), "solve with PINV fallback failed\n" );
    for ( integer i{0}; i < 3; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "fallback solve", i, 0, 5e-9 );

    base.Mv( x.data(), resid.data() );
    for ( integer i{0}; i < 3; ++i )
      check_close( resid[static_cast<size_t>(i)], rhs[static_cast<size_t>(i)], "fallback residual", i, 0, 5e-9 );
  }

} // namespace

int
main() {
  try {
    test_api_and_selection();
    test_zero_and_load_methods();
    test_sparse_mv_dup_and_print();
    test_check_matrix_and_tsolve();
    test_solver_combinations();
    test_last_solver_fallback();
    msg.green( "All extended BorderedCR coverage checks passed!\n" );
  } catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
    return 1;
  } catch ( ... ) {
    msg.red( "Unknown error\n" );
    return 1;
  }

  return 0;
}
