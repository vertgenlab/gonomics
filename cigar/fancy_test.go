package cigar

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
	"testing"
)

func TestMakeExplicit(t *testing.T) {
	// Sam-Formatted Cigar  5M2I1D5M
	// Reference Sequenced: ACGTA  CGTACG
	// Read Sequence:       ACTTATT GTACG
	// Explicit Cigar:      2M1Xt2M2Itt1D5M
	var expectedCigar string = "2M1Xt2M2Itt1D5M"

	var refSeq []dna.Base = dna.StringToBases("ACGTACGTACG")
	var readSeq []dna.Base = dna.StringToBases("ACTTATTGTACG")
	var c1 Cigar = Cigar{RunLength: 5, Op: 'M'}
	var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
	var c3 Cigar = Cigar{RunLength: 1, Op: 'D'}
	var c4 Cigar = Cigar{RunLength: 5, Op: 'M'}
	var samCigar []*Cigar = []*Cigar{&c1, &c2, &c3, &c4}

	explicitCigar := MakeExplicit(samCigar, readSeq, refSeq)

	fmt.Println("Input Cigar:", ToString(samCigar))
	fmt.Println("Output Cigar:", ToString(explicitCigar))

	reconstructedSeq := dna.BasesToString(GetExplicitSequence(explicitCigar, refSeq))

	fmt.Println("Original Sequence:", dna.BasesToString(readSeq))
	fmt.Println("Reconstructed Sequence:", reconstructedSeq)

	if strings.Compare(ToString(explicitCigar), expectedCigar) != 0 {
		log.Fatalln("Problem with conversion to explicit cigar")
	}
	if strings.Compare(dna.BasesToString(readSeq), reconstructedSeq) != 0 {
		log.Fatalln("Problem with sequence reconstruction from explicit cigar")
	}
}
