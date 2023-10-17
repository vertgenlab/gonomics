package pDna

// Float32Base encodes a DNA base as a probability vector using float32 precision.
// Note that gap probabilities can be stored implicitly as 1 - (A + C + G + T).
type Float32Base struct {
	A float32
	C float32
	G float32
	T float32
}
