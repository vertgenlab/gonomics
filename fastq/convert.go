package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
)

// ReverseComplement will modify the provided Fastq record
// to be the reverse complement of the input.  The name is
// unchanged, the sequences is reverse complemented, and the
// quality is reversed.
func ReverseComplement(record Fastq) {
	dna.ReverseComplement(record.Seq)
	ReverseQualUint8Record(record.Qual)
}
