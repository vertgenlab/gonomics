package pFasta

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/fasta"
	"strings"
	"log"
)

// AllAreEqual returns true if two []pFasta data structures are determined
// to be equal, false otherwise.
// Two []pFasta data structures are said to be equal if: (1) they are of the same
// length. (2) Each component pFasta struct in 'a' is equal to the pFasta struct
// at the corresponding index in 'b'. See the documentation for 'Equal' for a
// complete definition of the criterion for pFasta struct equality.
// Note that pFasta equality is parameterized by a user-specified relative 'precision'.
func AllAreEqual(a []PFasta, b []PFasta, precision float32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !IsEqual(a[i], b[i], precision) {
			return false
		}
	}
	return true
}

// IsEqual returns true if two input pFasta structs are 'equal', false otherwise.
// Two input pFasta structs are said to be equal if they have the same name,
// have the same number of bases in their sequence, and for each base, the
// base probabilities are equal within a user-specified level of relative precision.
func IsEqual(a PFasta, b PFasta, precision float32) bool {
	if strings.Compare(a.Name, b.Name) != 0 {
		return false
	}
	if len(a.Seq) != len(b.Seq) {
		return false
	}
	for i := range a.Seq {
		if !pDna.EqualBase(a.Seq[i], b.Seq[i], precision) {
			return false
		}
	}
	return true
}

// DistTrack reports a wig track from two input pFastas, providing base-by-base
// information about the similarity of the pFastas. Assumes the pFastas are aligned.
func DistTrack(a PFasta, b PFasta, outName string, defaultValue float64) wig.Wig {
	if len(a.Seq) != len(b.Seq) {
		log.Fatalf("Error (DistTrack): Input pFa sequences are not of the same length. len(a): %v. len(b): %v.", len(a.Seq), len(b.Seq))
	}

	if outName == "" {
		outName = a.Name
	}
	outWig := wig.Wig{StepType: "fixedStep", Chrom: outName, Start: 1, Step: 1, Span: 1, DefaultValue: defaultValue}
	outWig.Values = make([]float64, len(a.Seq))

	for pos := range len(a.Seq) {
		outWig.Values[pos] = pDna.Dist(a.Seq[pos], b.Seq[pos])
	}
	return outWig
}

// DistTrackFasta reports a wig track from an input pFasta and Fasta, providing base-by-base
// information about the similarity of the pFastas. Assumes the pFastas are aligned.
func DistTrackFasta(a PFasta, b fasta.Fasta, outName string, defaultValue float64) wig.Wig {
	return DistTrack(a, FaToPfa(b, 0, -1), outName, defaultValue)
}

// // TODO unimplemented and not used
// // DistTrackMulti reports a wig track from an input pFasta and Fasta, providing base-by-base
// // information about the similarity of the pFastas. Assumes the pFastas are aligned.
// func DistTrackMulti(a []PFasta, aName string, b []PFasta, bName string, outName string, defaultValue float64) wig.Wig {
// // 	a_len := 0
// 	for _, record := range a {
// 		if record.Name == aName {
// 			a_len = len(record.Seq)
// 		}
// 	}

// 	b_len := 0
// 	for _, record := range b {
// 		if record.Name == bName {
// 			b_len = len(record.Seq)
// 		}
// 	}

// 	a_single := Extract(a, 0, a_len, aName, aName, false)
// 	b_single := Extract(b, 0, b_len, aName, aName, false)
// 	return DistTrack(a_single, b_single, outName, defaultValue)

// }

