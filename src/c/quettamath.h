/*
 * QUETTAMATH - MATHEMATICAL LIBRARY FOR THE C PROGRAMMING LANGUAGE
 * COPYRIGHT (c) 2025 CYRIL JOHN MAGAYAGA
 */

#ifndef QUETTAMATH_H
#define QUETTAMATH_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double PI;
    double E;
    double TAU;
} qm_constant_t;
extern const qm_constant_t qm_constants;

double qm_positive(double x);
double qm_negative(double x);
double qm_add(size_t count, ...);
double qm_subtract(size_t count, ...);
double qm_multiply(size_t count, ...);
double qm_divide(size_t count, ...);
double qm_min(size_t count, ...);
double qm_max(size_t count, ...);
double qm_sq(double x);
double qm_sqrt(double x);
double qm_cb(double x);
double qm_cbrt(double x);
double qm_fabs(double x);
double qm_floor(double x);
double qm_ceil(double x);
double qm_trunc(double x);
long long qm_gcd(size_t count, ...);
long long qm_lcm(size_t count, ...);
double qm_degrees(double radians);
double qm_radians(double degrees);
static double qm_norm_angle(double x);
double qm_sin(double x);
double qm_cos(double x);
double qm_tan(double x);
double qm_csc(double x);
double qm_sec(double x);
double qm_cot(double x);

#ifndef positive
#define positive(x) qm_positive(x)
#endif

#ifndef gcd
#define gcd(count, ...) qm_gcd(count, __VA_ARGS__)
#endif

#ifdef __cplusplus
}
#endif

#endif // QUETTAMATH_H