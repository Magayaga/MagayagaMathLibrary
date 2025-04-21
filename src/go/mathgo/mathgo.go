/*
    MATHGO - MATHEMATICAL LIBRARY FOR THE GO PROGRAMMING LANGUAGE
	IT IS THE SIMILAR TO MATHMOJO

	Copyright (c) 2025 Cyril John Magayaga.
	It is originally written in Go programming language
*/
package mathgo

import (
	"fmt"
)

// Constants
type constants struct {
	// Mathematical constants
	PI float64
	E float64
	TAU float64
	goldenRatio float64
	PHI float64
	silverRatio float64

	// Physical constants
	C int64
	speedOfLight int64
	G float64
	gravitationalConstant float64
	e float64
	elementaryCharge float64
}

func constants() *constants {
	// PI
	pi_number := 3.141592653589793

	// EULER'S NUMBER
	e_number := 2.718281828459045

	// TAU IS EQUAL TO 2 * PI
	tau_number := 2 * pi_number

	// GOLDEN RATIO
	goldenRatio_number := 1.618033988749895
	phi_number := goldenRatio_number

	// SILVER RATIO
	silverRatio_number := 2.414213562373095

	// SPEED OF LIGHT
	speedOfLight_number := 299792458

	// Gravitational Constant
	gravitationalConstant_number := 0.0000000000667430

	// ELEMENTARY CHANGE
	elementaryCharge_number := 0.0000000000000000001602176634

	return &Constants{
		PI:  pi_number,
		E:   e_number,
		TAU: tau_number,
		goldenRatio: goldenRatio_number,
		PHI: phi_number,
		silverRatio: silverRatio_number,
		C: speedOfLight_number,
		speedOfLight: speedOfLight_number,
		G: gravitationalConstant_number
		gravitationalConstant: gravitationalConstant_number,
		e: elementaryCharge_number,
		elementaryCharge: elementaryCharge_number,
	}
}

// POSITIVE
func positive(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// NEGATIVE
func negative(x float64) float64 {
	if x > 0 {
		return -x
	}
	return x
}

// ADDITION
func add(values ...float64) float64 {
	sum := 0.0
	for _, v := range values {
		sum += v
	}
	return sum
}

// SUBTRACTION
func subtract(values ...float64) float64 {
	if len(values) == 0 {
		return 0.0
	}
	result := values[0]
	for _, v := range values[1:] {
		result -= v
	}
	return result
}

// MULTIPLICATION
func multiply(values ...float64) float64 {
	if len(values) == 0 {
		return 1.0
	}
	product := 1.0
	for _, v := range values {
		product *= v
	}
	return product
}

// DIVISION
func divide(values ...float64) float64 {
	if len(values) == 0 {
		return 0.0
	}
	result := values[0]
	for _, v := range values[1:] {
		if v == 0 {
			panic("Error: Division by zero")
		}
		result /= v
	}
	return result
}