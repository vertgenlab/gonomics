package vcf

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

var AncestorAlleles []dna.Base = []dna.Base{dna.A, dna.T, dna.G, dna.C, dna.A, dna.A, dna.C}

func TestVcfQueryAncestor(t *testing.T) {
	reader, _ := GoReadToChan("testdata/AncestorTest.vcf")
	var input []dna.Base
	var i int = 0
	for v := range reader {
		input = QueryAncestor(v)
		if input[0] != AncestorAlleles[i] {
			t.Errorf("Error in TestVcfQueryAncestor. Input: %s. Expected: %s.", dna.BaseToString(input[0]), dna.BaseToString(AncestorAlleles[i]))
		}
		i++
	}
}

func TestVcfAppendAncestor(t *testing.T) {
	reader, _ := GoReadToChan("testdata/Ancestor_No_Annotation.vcf")
	var i int = 0
	var input []dna.Base
	var allele []dna.Base = make([]dna.Base, 1)

	for v := range reader {
		allele[0] = AncestorAlleles[i]
		v = AppendAncestor(v, allele)
		input = QueryAncestor(v)
		if input[0] != AncestorAlleles[i] {
			t.Errorf("Error in TestVcfAppendAncestor. Input: %s. Expected: %s.", dna.BaseToString(input[0]), dna.BaseToString(AncestorAlleles[i]))
		}
		i++
	}
}

// This test just checks whether or not a line has been added. The print statements commented out allowed me to see that it was written in the right location.
// Although it seems like the header line order is not especially standardized, so I don't know if that part is important.
func TestAncestorFlagToHeader(t *testing.T) {
	_, header := GoReadToChan("testdata/Ancestor_No_Annotation.vcf")
	before := len(header.Text)
	header = AncestorFlagToHeader(header)
	after := len(header.Text)
	//DEBUG: PrintHeader(header)
	if after-before != 1 {
		t.Errorf("Error in TestAncestorFlagToHeader.")
	}

	//now we test a file that has no info columns in the header
	_, header = GoReadToChan("testdata/Ancestor_No_Info.vcf")
	before = len(header.Text)
	header = AncestorFlagToHeader(header)
	after = len(header.Text)
	//DEBUG: PrintHeader(header)
	if after-before != 1 {
		t.Errorf("Error in TestAncestorFlagToHeader.")
	}
}

var answers [][]dna.Base = [][]dna.Base{dna.StringToBases("A"), dna.StringToBases("T"), dna.StringToBases("A"), dna.StringToBases("CCT"), dna.StringToBases("N")}

func TestVcfAnnotateAncestorFromFa(t *testing.T) {
	reader, _ := GoReadToChan("testdata/TestVcfAnnotateAncestorFromFa.vcf")
	records := fasta.Read("testdata/testAncestorSequence.fa")
	var i int = 0
	var currRef, currAln int = 0, 0
	var currVcf Vcf

	for v := range reader {
		currVcf, currRef, currAln = AnnotateAncestorFromMultiFa(v, records, currRef, currAln)
		//DEBUG: fmt.Printf("Answer: %s. Expected:%s. RefPos: %v. AlnPos:%v.\n", dna.BasesToString(QueryAncestor(v)), dna.BasesToString(answers[i]), currRef, currAln)
		if dna.CompareSeqsIgnoreCase(QueryAncestor(currVcf), answers[i]) != 0 {
			t.Errorf("Error in TestVcfAnnotateAncestorFromFa. Expected: %s. Found: %s.", dna.BasesToString(answers[i]), dna.BasesToString(QueryAncestor(v)))
		}
		i++
	}
}

var IsRefAnswers []bool = []bool{true, false, false}
var IsAltAnswers []bool = []bool{false, true, false}

func IsAncestorTests(t *testing.T) {
	records, _ := Read("testdata/IsAncestor.vcf")
	for i, v := range records {
		if IsRefAncestor(v) != IsRefAnswers[i] {
			t.Errorf("Error in IsRefAncestor. Expected: %t. Found: %t.", IsRefAnswers[i], IsRefAncestor(v))
		}
		if IsAltAncestor(v) != IsAltAnswers[i] {
			t.Errorf("Error in IsAltAncestor. Expected: %t. Found: %t.", IsAltAnswers[i], IsAltAncestor(v))
		}
	}
}
