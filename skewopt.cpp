#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <string.h>

struct PolyScore {
  double skew;
  double alpha;
  double size_score;
  double combined_score;
  unsigned int num_real_roots;
  int valid;
};

// msieve_poly.cpp
void analyze_one_poly_hook( long, mpz_t*, mpz_t*, double, double*, double*, double*, unsigned int* );

void InitCoeffs( mpz_t*, mpz_t* );
void ClearCoeffs( mpz_t*, mpz_t* );
void msieve_compute_params( mpz_t*, mpz_t*, long, double*, double*, double*, double*, long* );
long deg( mpz_t*, long );
PolyScore ScoreAtSkew( long, mpz_t*, mpz_t*, double );
PolyScore OptimizeSkew( long, mpz_t*, mpz_t*, double );
double NominalSkew( long, mpz_t* );

const int MAX_ACOEFF = 10;
const double MIN_OPTIMIZED_SKEW = 1e-12;
const double MAX_OPTIMIZED_SKEW = 1e12;
const double SKEW_EXPANSION_FACTOR = 10.0;
const double LOG_SKEW_TOLERANCE = 1e-8;
const double INVERSE_GOLDEN_RATIO = 0.6180339887498948482;

int main(int argc, char* argv[] )
{
  int degree9_mode = ( argc > 1 && strcmp( argv[1], "-deg9" ) == 0 );
  int first_value = degree9_mode ? 2 : 1;
  int max_degree = degree9_mode ? 9 : 8;
  int num_coeffs = max_degree + 1;
  int argc_without_skew = first_value + 2 + num_coeffs;

  if ( argc != argc_without_skew && argc != argc_without_skew + 1 ) {
    printf("\n");
    printf("skewopt Version 1.0\n\n");
    printf("Usage: \"skewopt y0 y1 c0 c1 c2 c3 c4 c5 c6 c7 c8 [skew]\"");
    printf("\n");
    printf("       \"skewopt -deg9 y0 y1 c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 [skew]\"");
    printf("\n\n");
    return -1;
  }

  double requested_skew = 0.0;
  if ( argc == argc_without_skew + 1 ) {
    char* end = NULL;
    requested_skew = strtod( argv[argc_without_skew], &end );
    if ( end == argv[argc_without_skew] || *end != '\0' || !isfinite( requested_skew ) || requested_skew <= 0.0 ) {
      fprintf( stderr, "Invalid skew '%s': skew must be a positive finite number.\n", argv[argc_without_skew] );
      return -1;
    }
  }

  mpz_t c[MAX_ACOEFF];
  mpz_t y[2];

  InitCoeffs( c, y );

  int valid_input = mpz_set_str( y[0], argv[first_value], 10 ) == 0 &&
                    mpz_set_str( y[1], argv[first_value + 1], 10 ) == 0;

  for ( int i = 0; i < num_coeffs; i++ )
    valid_input = mpz_set_str( c[i], argv[first_value + 2 + i], 10 ) == 0 &&
                  valid_input;

  if ( !valid_input ) {
    fprintf( stderr, "Invalid polynomial: all coefficients must be base-10 integers.\n" );
    ClearCoeffs( c, y );
    return -1;
  }

  double skew = 0.0;
  double alpha = 0.0;
  double size_score = 0.0;
  double combined_score = 0.0;
  long num_real_roots = 0;

  long degree = deg( c, max_degree );
  if ( degree < 1 ) {
    fprintf( stderr, "Invalid polynomial: the algebraic polynomial must have degree at least 1.\n" );
    ClearCoeffs( c, y );
    return -1;
  }
  if ( argc == argc_without_skew + 1 ) {
    unsigned int fixed_num_real_roots = 0;
    skew = requested_skew;
    analyze_one_poly_hook( degree, c, y, skew, &size_score, &alpha, &combined_score, &fixed_num_real_roots );
    num_real_roots = fixed_num_real_roots;
    printf("Skew: %.12g\n", skew );
  }
  else {
    msieve_compute_params( c, y, degree, &skew, &alpha, &size_score, &combined_score, &num_real_roots );
    printf("Best Skew: %.9g\n", skew );
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

  double nominal_skew = NominalSkew( degree, c );
  if ( !isfinite( nominal_skew ) || nominal_skew <= 0.0 )
    nominal_skew = 1.0;
  PolyScore best = OptimizeSkew( degree, c, y, nominal_skew );

  *skew = best.skew;
  *size_score = best.size_score;
  *alpha = best.alpha;
  *combined_score = best.combined_score;
  *num_real_roots = best.num_real_roots;
}

PolyScore ScoreAtSkew( long degree, mpz_t* c, mpz_t* y, double skew ) {
  PolyScore result = { skew, 0.0, 0.0, -INFINITY, 0, 0 };
  analyze_one_poly_hook( degree, c, y, skew, &result.size_score,
                         &result.alpha, &result.combined_score,
                         &result.num_real_roots );
  result.valid = isfinite( result.combined_score );
  return result;
}

static void KeepBetterScore( PolyScore* best, const PolyScore* candidate ) {
  if ( candidate->valid &&
       ( !best->valid || candidate->combined_score > best->combined_score ) )
    *best = *candidate;
}

static double ExpandedUpperSkew( double skew, double upper_limit ) {
  if ( skew >= upper_limit / SKEW_EXPANSION_FACTOR )
    return upper_limit;
  return skew * SKEW_EXPANSION_FACTOR;
}

// Expand by decades until the maximum is bracketed, then refine in log-skew
// space. Murphy E is assumed to have one maximum, as in the original search.
// If the nominal skew lies outside the search limits, score it separately so
// optimization can never make its score worse.
PolyScore OptimizeSkew( long degree, mpz_t* c, mpz_t* y, double nominal_skew ) {
  double lower_limit = MIN_OPTIMIZED_SKEW;
  double upper_limit = MAX_OPTIMIZED_SKEW;
  PolyScore nominal = ScoreAtSkew( degree, c, y, nominal_skew );
  PolyScore best = nominal;
  double center_skew = fmin( upper_limit,
                             fmax( lower_limit, nominal_skew ) );
  PolyScore center = center_skew == nominal_skew
                       ? nominal : ScoreAtSkew( degree, c, y, center_skew );
  KeepBetterScore( &best, &center );

  double left_skew = fmax( lower_limit,
                           center_skew / SKEW_EXPANSION_FACTOR );
  double right_skew = ExpandedUpperSkew( center_skew, upper_limit );
  PolyScore left = left_skew == center_skew
                     ? center : ScoreAtSkew( degree, c, y, left_skew );
  PolyScore right = right_skew == center_skew
                      ? center : ScoreAtSkew( degree, c, y, right_skew );
  KeepBetterScore( &best, &left );
  KeepBetterScore( &best, &right );

  while ( true ) {
    int improve_left = left.skew < center.skew && left.valid &&
                       ( !center.valid ||
                         left.combined_score > center.combined_score );
    int improve_right = right.skew > center.skew && right.valid &&
                        ( !center.valid ||
                          right.combined_score > center.combined_score );
    if ( !improve_left && !improve_right )
      break;

    if ( improve_left &&
         ( !improve_right ||
           left.combined_score >= right.combined_score ) ) {
      if ( left.skew == lower_limit ) {
        right = center;
        center = left;
        break;
      }
      right = center;
      center = left;
      left_skew = fmax( lower_limit,
                        center.skew / SKEW_EXPANSION_FACTOR );
      left = ScoreAtSkew( degree, c, y, left_skew );
      KeepBetterScore( &best, &left );
    }
    else {
      if ( right.skew == upper_limit ) {
        left = center;
        center = right;
        break;
      }
      left = center;
      center = right;
      right_skew = ExpandedUpperSkew( center.skew, upper_limit );
      right = ScoreAtSkew( degree, c, y, right_skew );
      KeepBetterScore( &best, &right );
    }
  }

  double log_left = log( left.skew );
  double log_right = log( right.skew );
  double log_first = log_right - INVERSE_GOLDEN_RATIO *
                     ( log_right - log_left );
  double log_second = log_left + INVERSE_GOLDEN_RATIO *
                      ( log_right - log_left );
  PolyScore first = ScoreAtSkew( degree, c, y, exp( log_first ) );
  PolyScore second = ScoreAtSkew( degree, c, y, exp( log_second ) );
  KeepBetterScore( &best, &first );
  KeepBetterScore( &best, &second );

  while ( log_right - log_left > LOG_SKEW_TOLERANCE ) {
    if ( first.combined_score < second.combined_score ) {
      log_left = log_first;
      log_first = log_second;
      first = second;
      log_second = log_left + INVERSE_GOLDEN_RATIO *
                   ( log_right - log_left );
      second = ScoreAtSkew( degree, c, y, exp( log_second ) );
      KeepBetterScore( &best, &second );
    }
    else {
      log_right = log_second;
      log_second = log_first;
      second = first;
      log_first = log_right - INVERSE_GOLDEN_RATIO *
                  ( log_right - log_left );
      first = ScoreAtSkew( degree, c, y, exp( log_first ) );
      KeepBetterScore( &best, &first );
    }
  }
  return best;
}

// Compute the nominal skew
double NominalSkew( long degree, mpz_t* c ) {
  return pow( fabs( mpz_get_d( c[0] ) / mpz_get_d( c[degree] ) ), 1.0 / (double)degree );
}
