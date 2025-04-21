/*
 * QUETTAMATH - MATHEMATICAL LIBRARY FOR THE C PROGRAMMING LANGUAGE
 * COPYRIGHT (c) 2025 CYRIL JOHN MAGAYAGA
 */

#ifndef QUETTAMATH_H
#define QUETTAMATH_H

#include <stddef.h>

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
double qm_fabs(double x);
double qm_floor(double x);
double qm_ceil(double x);
double qm_trunc(double x);
long long qm_gcd(size_t count, ...);
long long qm_lcm(size_t count, ...);
double qm_degrees(double radians);
double qm_radians(double degrees);

#endif // QUETTAMATH_H