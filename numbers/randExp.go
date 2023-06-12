package numbers

import (
	"math"
	"math/rand"
)

// RandGeometric returns a geometrically-distributed random variate using the inverse transform.
// Note the geometric distribution has CDF: F(x) = 1 - (1-p)^x.
// The inverse transform is derived as (1-p)^x = 1 - F(x) -> x*log(1-p) = log(1 - F(x)) -> x = log(1 - F(x)) / log(1-p).
// Note that this is the version of the geometric distribution with support from 0 to +INF.
func RandGeometric(p float64) int {
	r := rand.Float64()
	return int(math.Floor(math.Log(1-r) / math.Log(1-p)))
}

// RandExp Returns a random variable as a float64 from a standard exponential distribution. f(x)=e**-x.
// Algorithm from Ahrens, J.H. and Dieter, U. (1972). Computer methods for sampling from the exponential and normal distributions. Comm. ACM, 15, 873-882.
func RandExp() (float64, float64) {
	//q series where q[k-1] = sum(log(2)^k / k!) for k=1,2,...n
	q := [16]float64{0.6931471805599453, 0.9333736875190459, 0.9888777961838675, 0.9984959252914960, 0.9998292811061389, 0.9999833164100727, 0.9999985691438767, 0.9999998906925558, 0.9999999924734159, 0.9999999995283275, 0.9999999999728814, 0.9999999999985598, 0.9999999999999289, 0.9999999999999968, 0.9999999999999999, 1.0000000000000000}

	var a float64 = 0.0
	var r float64 = rand.Float64()
	for r <= 0.0 || r >= 1.0 { //prevents the case where u is exactly 0 or 1, which breaks the code.
		r = rand.Float64()
	}

	for 1 > 0 {
		r += r
		if r > 1.0 {
			break
		}
		a += q[0]
	}
	r -= 1
	if r <= q[0] {
		return a + r, ExpDist(a + r)
	}

	var i int = 0
	ustart := rand.Float64()
	umin := ustart

	for r > q[i] {
		ustart = rand.Float64()
		if umin > ustart {
			umin = ustart
		}
		i++
	}
	return a + umin*q[0], ExpDist(a + umin*q[0])
}
