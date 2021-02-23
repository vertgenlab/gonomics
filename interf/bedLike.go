package interf

// bedLike interface compatible with any data type that
// implements methods to retrieve Chr, start and end pos
// Current data types that implement these methods:
// bed.Bed
// axt.Axt
// vcf.Vcf

type BedLike interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
}
