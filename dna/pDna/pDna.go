package pDna

import (
	"log"
	"math"
)

// Float32Base encodes a DNA base as a probability vector using float32 precision.
// gap probabilities can be stored implicitly as 1 - (A + C + G + T).
type Float32Base struct {
	A float32
	C float32
	G float32
	T float32
}

// Uint8Base encodes a DNA base as a probability vector using uint8 probability encoding.
// A uint8 value x can be converted to a probability as follows: p = 1 / (2^(x/10)).
// gap probabilities can be stored implicitly as 1 - (A + C + G + T).
type Uint8Base struct {
	A uint8
	C uint8
	G uint8
	T uint8
}

// Float32BaseToUint8Base takes an input Float32Base and converts it to a Uint8Base, where
// each base probability is represented by an encoded uint8.
func Float32BaseToUint8Base(base Float32Base) Uint8Base {
	return Uint8Base{
		A: Float32ProbToUint8(base.A),
		C: Float32ProbToUint8(base.C),
		G: Float32ProbToUint8(base.G),
		T: Float32ProbToUint8(base.T),
	}
}

// Float32ProbToUint8 converts a probability in float32 to an encoded probability in uint8.
func Float32ProbToUint8(prob float32) uint8 {
	if prob < 0 || prob > 1 {
		log.Fatalf("Error: float probability must be between 0 and 1. Found: %v.\n", prob)
	}
	if prob == 0 {
		return 0
	}
	return uint8(math.Round(-10 * math.Log2(float64(prob))))
}

// Uint8ToFloat32Prob converts an encoded probability in uint8 format into a probability in float32 format.
func Uint8ToFloat32Prob(u uint8) float32 {
	if u == 0 {
		return 0
	}
	return float32(1 / math.Pow(2, float64(u)/10.0))
}

// MakeUint8ToFloat32Map makes a map relating uint8 probabilities to the corresponding float32
// representations
func MakeUint8ToFloat32Map() map[uint8]float32 {
	var answer = make(map[uint8]float32)
	var u uint8
	for u = 0; u < 255; u++ {
		answer[u] = Uint8ToFloat32Prob(u)
	}
	answer[255] = Uint8ToFloat32Prob(255)
	return answer
}
