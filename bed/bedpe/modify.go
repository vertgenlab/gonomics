package bedpe

import "github.com/vertgenlab/gonomics/bed"

//contactToMidpoint returns a bedpe where the contact has been collapsed to its midpoint
func contactToMidpoint(bp BedPe) BedPe {
	aMidpoint := (bp.A.ChromStart + bp.A.ChromEnd) / 2
	bMidpoint := (bp.B.ChromStart + bp.B.ChromEnd) / 2
	A := bed.Bed{Chrom: bp.A.Chrom, ChromStart: aMidpoint, ChromEnd: aMidpoint + 1, FieldsInitialized: bp.A.FieldsInitialized}
	B := bed.Bed{Chrom: bp.B.Chrom, ChromStart: bMidpoint, ChromEnd: bMidpoint + 1, FieldsInitialized: bp.A.FieldsInitialized}
	return BedPe{A, B}
}

//contactsToMidpoints returns a bedpe where the contacts have been collapsed to their midpoints
func contactsToMidpoints(bps []BedPe) []BedPe {
	var answers = make([]BedPe, len(bps))
	for i := range bps {
		answers[i] = contactToMidpoint(bps[i])
	}
	return answers
}
