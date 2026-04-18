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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#include <cmath>

static Utils::Console msg( &std::cout, 4 );

using real_type = double;
using integer   = alglin::integer;
using BCR       = alglin::BorderedCR<real_type>;
using MatW      = alglin::MatrixWrapper<real_type>;
using Matrix    = alglin::Matrix<real_type>;

static
void
check_close(
  real_type value,
  real_type expected,
  char const where[],
  integer i,
  integer j
) {
  UTILS_ASSERT(
    std::abs(value-expected) < 1e-12,
    "{} mismatch at ({},{}): found {}, expected {}\n",
    where, i, j, value, expected
  );
}

static
void
test_matw_views() {
  Utils::ThreadPool1 TP(4);
  BCR                bcr(&TP);

  integer constexpr nblock{2};
  integer constexpr n     {2};
  integer constexpr qr    {1};
  integer constexpr qx    {1};
  integer constexpr nr    {2};
  integer constexpr nx    {1};

  bcr.allocate( nblock, n, qr, qx, nr, nx );
  bcr.fill_zero();

  Matrix bottom_storage( nr+2, bcr.Nc() );
  MatW   bottom_view;
  bottom_storage.fill( -999 );
  bottom_storage.view_block( 1, 0, nr, bcr.Nc(), bottom_view );

  for ( integer j{0}; j < bottom_view.ncols(); ++j )
    for ( integer i{0}; i < bottom_view.nrows(); ++i )
      bottom_view(i,j) = 100 + 10*j + i;

  bcr.load_bottom2( bottom_view );

  for ( integer j{0}; j < n; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.C(0,i,j), bottom_view(i,j), "load_bottom2/C0", i, j );

  for ( integer j{0}; j < n; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.C(nblock,i,j), bottom_view(i,n+j), "load_bottom2/CN", i, j );

  for ( integer j{0}; j < qx; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.Cq(i,j), bottom_view(i,2*n+j), "load_bottom2/Cq", i, j );

  for ( integer j{0}; j < nx; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.F(i,j), bottom_view(i,2*n+qx+j), "load_bottom2/F", i, j );

  bcr.fill_zero();

  Matrix c2f_storage( nr+1, 2*n+nx );
  MatW   c2f_view;
  c2f_storage.fill( -777 );
  c2f_storage.view_block( 0, 0, nr, 2*n+nx, c2f_view );

  for ( integer j{0}; j < c2f_view.ncols(); ++j )
    for ( integer i{0}; i < c2f_view.nrows(); ++i )
      c2f_view(i,j) = 200 + 10*j + i;

  bcr.add_to_C2F( 1, c2f_view );

  for ( integer j{0}; j < n; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.C(1,i,j), c2f_view(i,j), "add_to_C2F/C1", i, j );

  for ( integer j{0}; j < n; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.C(2,i,j), c2f_view(i,n+j), "add_to_C2F/C2", i, j );

  for ( integer j{0}; j < nx; ++j )
    for ( integer i{0}; i < nr; ++i )
      check_close( bcr.F(i,j), c2f_view(i,2*n+j), "add_to_C2F/F", i, j );
}

static
void
test_parallel_reconfigure_and_info() {
  Utils::ThreadPool1 TP(4);
  BCR                bcr(&TP,1);

  bcr.allocate( 16, 2, 0, 0, 0, 0 );

  std::string info{ bcr.info() };
  UTILS_ASSERT0( info.find("nblk:1") != std::string::npos, "expected nblk:1 in initial info\n" );

  bcr.set_num_parallel_block(4);
  info = bcr.info();
  UTILS_ASSERT0( info.find("nblk:4") != std::string::npos, "expected nblk:4 after set_num_parallel_block(4)\n" );

  bcr.set_num_parallel_block(2);
  info = bcr.info();
  UTILS_ASSERT0( info.find("nblk:2") != std::string::npos, "expected nblk:2 after set_num_parallel_block(2)\n" );

  bcr.select_last_LU();
  bcr.select_last2_PINV();
  info = bcr.info();
  UTILS_ASSERT0(
    info.find("last   = LastBlock LU") != std::string::npos,
    "expected primary last-block solver in info()\n"
  );
  UTILS_ASSERT0(
    info.find("last2  = LastBlock PINV") != std::string::npos,
    "expected fallback last-block solver in info()\n"
  );
}

int
main() {
  try {
    test_matw_views();
    test_parallel_reconfigure_and_info();
    msg.green( "All BorderedCR regression checks passed!\n" );
  } catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
    return 1;
  } catch ( ... ) {
    msg.red( "Unknown error\n" );
    return 1;
  }

  return 0;
}
