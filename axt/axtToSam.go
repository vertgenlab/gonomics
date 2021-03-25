package axt

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

// ToSam converts an Axt record into a Sam (sam.SamAln) record
// TODO: Add logic to add hard clip on ends of axt that contain no sequence alignment
func ToSam(axtFmt Axt) *sam.SamAln {
	var answer *sam.SamAln = &sam.SamAln{
		QName: axtFmt.QName,
		Flag:  setStrandFlag(axtFmt.QStrandPos),
		RName: axtFmt.RName,
		Pos:   int64(axtFmt.RStart),
		MapQ:  255, // mapping quality setting to 255 because we are not calculating it
		Cigar: PairSeqToCigar(axtFmt.RSeq, axtFmt.QSeq),
		RNext: "*",
		PNext: 0,
		TLen:  int64(axtFmt.REnd - axtFmt.RStart),
		Seq:   dna.RemoveBase(axtFmt.QSeq, dna.Gap),
		Qual:  "*",
		Extra: fmt.Sprintf("AS:i:%d\tXS:i:%d\tXE:i:%d", axtFmt.Score, axtFmt.QStart, axtFmt.QEnd),
		//AS=alignment score, XS=query start , XE= query end position
	}
	//answer.Cigar = append(answer.Cigar, &cigar.Cigar{Op: 'H', RunLength: axtFmt.QStart-1} )
	return answer
}

// PairSeqToCigar takes in two sequences representing an alignment between the two
// sequences and outputs a slice of cigars describing the input alignment.
func PairSeqToCigar(a []dna.Base, b []dna.Base) []*cigar.Cigar {
	var align []*cigar.Cigar = make([]*cigar.Cigar, 0)
	curr := &cigar.Cigar{}
	var i int
	for i = 0; i < len(a); i++ {
		switch true {
		case a[i] != dna.Gap && b[i] != dna.Gap && a[i] == b[i]: //match, bases equal
			curr = equalMatchCigar(a, b, i)
			i += curr.RunLength - 1
			align = append(align, curr)
		case a[i] != dna.Gap && b[i] != dna.Gap && a[i] != b[i]:
			curr = diffMatchCigar(a, b, i)
			i += curr.RunLength - 1
			align = append(align, curr)
		case a[i] == dna.Gap && b[i] != dna.Gap: //target has gap, non-gap base in target
			curr = insertCigar(a, b, i)
			i += curr.RunLength - 1
			align = append(align, curr)
		case a[i] != dna.Gap && b[i] == dna.Gap: //query is a gap, contains sequence
			curr = deletionCigar(a, b, i)
			i += curr.RunLength - 1
			align = append(align, curr)
		default:
			log.Fatalf("Error: not catching case when bases are %c and %c\n", dna.BaseToRune(a[i]), dna.BaseToRune(b[i]))
		}
	}
	return align
}

// equalMatchCigar is a helper function for PairSeqToCigar
func equalMatchCigar(a []dna.Base, b []dna.Base, index int) *cigar.Cigar {
	match := &cigar.Cigar{Op: '=', RunLength: 1}
	var i int
	for i = index + 1; i < len(a); i++ {
		if a[i] == b[i] && a[i] != dna.Gap && b[i] != dna.Gap {
			match.RunLength++
		} else {
			return match
		}
	}
	return match
}

// diffMatchCigar is a helper function for PairSeqToCigar
func diffMatchCigar(a []dna.Base, b []dna.Base, index int) *cigar.Cigar {
	match := &cigar.Cigar{Op: 'X', RunLength: 1}
	var i int
	for i = index + 1; i < len(a); i++ {
		if a[i] != b[i] && a[i] != dna.Gap && b[i] != dna.Gap {
			match.RunLength++
		} else {
			return match
		}
	}
	return match
}

// insertCigar is a helper function for PairSeqToCigar
func insertCigar(a []dna.Base, b []dna.Base, index int) *cigar.Cigar {
	insertion := &cigar.Cigar{Op: 'I', RunLength: 1}
	var i int
	//starting loop at +1 since we already checked in the wrapper function above
	for i = index + 1; i < len(a); i++ {
		if a[i] == dna.Gap {
			insertion.RunLength++
		} else {
			return insertion
		}
	}
	return insertion
}

// deletionCigar is a helper function for PairSeqToCigar
func deletionCigar(a []dna.Base, b []dna.Base, index int) *cigar.Cigar {
	deletion := &cigar.Cigar{Op: 'D', RunLength: 1}
	var i int
	//starting loop at +1 since we already checked in the wrapper function above
	for i = index + 1; i < len(a); i++ {
		if b[i] == dna.Gap {
			deletion.RunLength++
		} else {
			return deletion
		}
	}
	return deletion
}

// setStrandFlag is a helper function for ToSam which sets
// the correct bit in the flag based on if we are aligning
// to the positive or negative strand.
func setStrandFlag(strand bool) int64 {
	if strand {
		return 0
	} else {
		return 16
	}
}
