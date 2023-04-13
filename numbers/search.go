package numbers

import (
	"math"
)

// GoldenSectionMaxSearch returns returns a local maximum for the input range a to b from an input function f with a user specified error epsilon.
// For unimodal functions, we are guaranteed to find the maximum if it is contained within the interval a to b. For polymodal functions, this function will converge on one of its modes contained in the interval a to b.
func GoldenSectionMaxSearch(f func(float64) float64, a float64, b float64, epsilon float64) float64 {
	a, b = Min(a, b), Max(a, b) //set a to be less than b.
	c := b + (a-b)/math.Phi
	d := a + (b-a)/math.Phi
	fc := f(c)
	fd := f(d)

	for math.Abs(b-a) > epsilon {
		if fc < fd { //reject left
			a = c
			c = d
			fc = fd
			d = a + (b-a)/math.Phi
			fd = f(d)
		} else {
			b = d
			d = c
			fd = fc
			c = b + (a-b)/math.Phi
			fc = f(c)
		}
	}
	return (a + b) / 2.0
}

// GoldenSectionMinSearch returns a local minimum for the input range a to b from an input function f with a user specified error epsilon.
// For a unimodal function, we are guaranteed to find the minimum if it is contained within the interval a to b. For polymodal functions, this function will converge on one of its modes contained in the interval a to b.
func GoldenSectionMinSearch(f func(float64) float64, a float64, b float64, epsilon float64) float64 {
	a, b = Min(a, b), Max(a, b) //set a to be less than b.
	c := b + (a-b)/math.Phi
	d := a + (b-a)/math.Phi
	fc := f(c)
	fd := f(d)
	for math.Abs(b-a) > epsilon {
		if fc < fd { //reject right
			b = d
			d = c
			fd = fc
			c = b + (a-b)/math.Phi
			fc = f(c)
		} else {
			a = c
			c = d
			fc = fd
			d = a + (b-a)/math.Phi
			fd = f(d)
		}
	}
	return (c + d) / 2.0
}
