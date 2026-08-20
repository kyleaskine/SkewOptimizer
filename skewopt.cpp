#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <string.h>

// msieve_poly.cpp
void analyze_one_poly_hook( long, mpz_t*, mpz_t*, double, double*, double*, double*, unsigned int* );

void InitCoeffs( mpz_t*, mpz_t* );
void ClearCoeffs( mpz_t*, mpz_t* );
void msieve_compute_params( mpz_t*, mpz_t*, long, double*, double*, double*, double*, long* );
long deg( mpz_t*, long );
double ComputeMaxima( long, mpz_t*, mpz_t*, double, double, double, double, double, double );
double NominalSkew( long, mpz_t* );

const int MAX_ACOEFF = 9;

int main(int argc, char* argv[] )
{
  if ( argc != 12 && argc != 13 ) {
    printf("\n");
    printf("skewopt Version 1.0\n\n");
    printf("Usage: \"skewopt y0 y1 c0 c1 c2 c3 c4 c5 c6 c7 c8 [skew]\"");
    printf("\n\n");
    return -1;
  }

  double requested_skew = 0.0;
  if ( argc == 13 ) {
    char* end = NULL;
    requested_skew = strtod( argv[12], &end );
    if ( end == argv[12] || *end != '\0' || !isfinite( requested_skew ) || requested_skew <= 0.0 ) {
      fprintf( stderr, "Invalid skew '%s': skew must be a positive finite number.\n", argv[12] );
      return -1;
    }
  }

  mpz_t c[MAX_ACOEFF];
  mpz_t y[2];

  InitCoeffs( c, y );

  mpz_set_str( y[0], argv[1], 10 );
  mpz_set_str( y[1], argv[2], 10 );

  mpz_set_str( c[0], argv[3], 10 );
  mpz_set_str( c[1], argv[4], 10 );
  mpz_set_str( c[2], argv[5], 10 );
  mpz_set_str( c[3], argv[6], 10 );
  mpz_set_str( c[4], argv[7], 10 );
  mpz_set_str( c[5], argv[8], 10 );
  mpz_set_str( c[6], argv[9], 10 );
  mpz_set_str( c[7], argv[10], 10 );
  mpz_set_str( c[8], argv[11], 10 );

  double skew = 0.0;
  double alpha = 0.0;
  double size_score = 0.0;
  double combined_score = 0.0;
  double anorm = 0.0;
  double rnorm = 0.0;
  long num_real_roots = 0;

  long degree = deg( c, MAX_ACOEFF - 1 );
  if ( argc == 13 ) {
    unsigned int fixed_num_real_roots = 0;
    skew = requested_skew;
    analyze_one_poly_hook( degree, c, y, skew, &size_score, &alpha, &combined_score, &fixed_num_real_roots );
    num_real_roots = fixed_num_real_roots;
    printf("Skew: %.12g\n", skew );
  }
  else {
    msieve_compute_params( c, y, degree, &skew, &alpha, &size_score, &combined_score, &num_real_roots );
    printf("Best Skew: %.5f\n", skew );
  }
  printf("MurphyE: %.8e\n", combined_score );

  ClearCoeffs( c, y );

  return 0;
}

void InitCoeffs( mpz_t* c, mpz_t* y ) {

  long i;
  for ( i = 0; i < MAX_ACOEFF; i++ )
    mpz_init( c[i] );

  for ( i = 0; i < 2; i++ )
    mpz_init( y[i] );
}

void ClearCoeffs( mpz_t* c, mpz_t* y ) {

  long i;
  for ( i = 0; i < MAX_ACOEFF; i++ )
    mpz_clear( c[i] );

  for ( i = 0; i < 2; i++ )
    mpz_clear( y[i] );
}


long deg( mpz_t* c, long maxdeg ) {

  long degree = maxdeg;
  for ( ; degree >= 0 ; degree-- ) {
    if ( mpz_cmp_ui( c[degree], 0 ) != 0 )
      break;
  }
  return degree;
}

// use msieve library
void msieve_compute_params( mpz_t* c, mpz_t* y, long degree, double* skew, double* alpha, double* size_score, double* combined_score, long* num_real_roots ) {

  double best_skewness = NominalSkew( degree, c );
  // Use a precision target relative to the skew magnitude.  An absolute
  // target (the old 0.0000001) is unreachable for large skews: near a skew
  // of ~1e9 one ULP (~2.4e-7) already exceeds it, so the bisection interval
  // can never shrink below it.  A relative target resolves the skew to ~9
  // significant figures at any magnitude.
  double min_precision = best_skewness * 1e-9;
  double tmpSkew = ComputeMaxima( degree, c, y, 0.1*best_skewness, -1.0, -1.0, 9.0*best_skewness, -1.0, min_precision );
  if ( tmpSkew > 0.0 )
    best_skewness = tmpSkew;

  double local_size_score = -1.0;
  double root_score = -1.0;
  double local_combined_score = -1.0;
  unsigned int local_num_real_roots = -1;
  analyze_one_poly_hook( degree, c, y, best_skewness, &local_size_score, &root_score, &local_combined_score, &local_num_real_roots );
  
  *skew = best_skewness;
  *size_score = local_size_score;
  *alpha = root_score;
  *combined_score = local_combined_score;
  *num_real_roots = local_num_real_roots;
}

// Find the maxima for combined_score (assumption: only one maxima)
double ComputeMaxima( long degree, mpz_t* c, mpz_t* y, double lowerbound, double lower_score, double middle_score, double upperbound, double upper_score, double min_precision ) {

  double size_score = 0.0;
  double root_score = 0.0;
  unsigned int num_real_roots = 0;

  double middle        = ( lowerbound + upperbound ) / 2.0;
  double bottomquarter = ( lowerbound + middle ) / 2.0;
  double topquarter    = ( middle + upperbound ) / 2.0;
  double precision     = ( upperbound - lowerbound ) / 4.0;

  // The interval can only be bisected further if its quarter points are all
  // distinct in double precision.  Once they collapse onto the endpoints the
  // bracket cannot narrow, so we must stop regardless of the requested
  // precision -- otherwise recursing on an unchanged interval loops forever.
  // (This is what hangs large skews: at ~1e9 the interval freezes 2 ULP wide,
  // where (upperbound-lowerbound)/2 rounds [bottomquarter,topquarter] right
  // back to [lowerbound,upperbound].)
  int subdividable = ( lowerbound < bottomquarter && bottomquarter < middle &&
                       middle < topquarter && topquarter < upperbound );

  if ( lower_score < 0.0 ) // uninitialized
    analyze_one_poly_hook( degree, c, y, lowerbound, &size_score, &root_score, &lower_score, &num_real_roots );

  if ( middle_score < 0.0 ) // uninitialized
    analyze_one_poly_hook( degree, c, y, middle, &size_score, &root_score, &middle_score, &num_real_roots );

  double bottomquarter_score = 0.0;
  analyze_one_poly_hook( degree, c, y, bottomquarter, &size_score, &root_score, &bottomquarter_score, &num_real_roots );

  if ( ( lower_score > bottomquarter_score ) || ( bottomquarter_score > middle_score ) ) {
    if ( precision < min_precision || !subdividable )
      return lower_score > bottomquarter_score ? lowerbound : ( bottomquarter_score > middle_score ? bottomquarter : middle );

    return ComputeMaxima( degree, c, y, lowerbound, lower_score, bottomquarter_score, middle, middle_score, min_precision );
  }

  if ( upper_score < 0.0 ) // uninitialized
    analyze_one_poly_hook( degree, c, y, upperbound, &size_score, &root_score, &upper_score, &num_real_roots );

  double topquarter_score = 0.0;
  analyze_one_poly_hook( degree, c, y, topquarter, &size_score, &root_score, &topquarter_score, &num_real_roots );

  if ( precision < min_precision || !subdividable )
    return middle_score > topquarter_score ? middle : ( topquarter_score > upper_score ? topquarter : upperbound );

  if ( middle_score > topquarter_score )
    return ComputeMaxima( degree, c, y, bottomquarter, bottomquarter_score, middle_score, topquarter, topquarter_score, min_precision );

  return ComputeMaxima( degree, c, y, middle, middle_score, topquarter_score, upperbound, upper_score, min_precision );
}

// Compute the nominal skew
double NominalSkew( long degree, mpz_t* c ) {
  return pow( fabs( mpz_get_d( c[0] ) / mpz_get_d( c[degree] ) ), 1.0 / (double)degree );
}
