package _interface

// bedLike interface compatible with any data type that
// implements methods to retrieve Chr, start and end pos
// Current data types that implement these methods:
// bed.Bed

type BedLike interface {
	GetChr() string
	GetStart() int
	GetEnd() int
}
