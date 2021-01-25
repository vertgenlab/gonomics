package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
	//DEBUG: "fmt"
)

var AncestorAlleles []dna.Base = []dna.Base{dna.A, dna.T, dna.G, dna.C, dna.A, dna.A, dna.C}

func TestVcfQueryAncestor(t *testing.T) {
	reader, _ := GoReadToChan("testdata/AncestorTest.vcf")
	var input []dna.Base
	var i int = 0
	for v := range reader {
		input = VcfQueryAncestor(v)
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
		VcfAppendAncestor(v, allele)
		input = VcfQueryAncestor(v)
		if input[0] != AncestorAlleles[i] {
			t.Errorf("Error in TestVcfAppendAncestor. Input: %s. Expected: %s.", dna.BaseToString(input[0]), dna.BaseToString(AncestorAlleles[i]))
		}
		i++
	}
}

//This test just checks whether or not a line has been added. The print statements commented out allowed me to see that it was written in the right location.
//Although it seems like the header line order is not especially standardized, so I don't know if that part is important.
func TestAncestorFlagToHeader(t *testing.T) {
	_, header := GoReadToChan("testdata/Ancestor_No_Annotation.vcf")
	before := len(header.Text)
	AncestorFlagToHeader(header)
	after := len(header.Text)
	//DEBUG: PrintHeader(header)
	if after-before != 1 {
		t.Errorf("Error in TestAncestorFlagToHeader.")
	}

	//now we test a file that has no info columns in the header
	_, header = GoReadToChan("testdata/Ancestor_No_Info.vcf")
	before = len(header.Text)
	AncestorFlagToHeader(header)
	after = len(header.Text)
	//DEBUG: PrintHeader(header)
	if after-before != 1 {
		t.Errorf("Error in TestAncestorFlagToHeader.")
	}
}

var answers [][]dna.Base = [][]dna.Base{{dna.A}, {dna.T}, {dna.A}, {dna.C, dna.C, dna.T}, {dna.N}}

func TestVcfAnnotateAncestorFromFa(t *testing.T) {
	reader, _ := GoReadToChan("testdata/TestVcfAnnotateAncestorFromFa.vcf")
	records := fasta.Read("testdata/testAncestorSequence.fa")
	var i int = 0

	for v := range reader {
		VcfAnnotateAncestorFromFa(v, records)
		//DEBUG: fmt.Printf("Answer: %s. Expected:%s. \n", dna.BasesToString(GVcfQueryAncestor(v)), dna.BasesToString(answers[i]))
		if dna.CompareSeqsIgnoreCase(VcfQueryAncestor(v), answers[i]) != 0 {
			t.Errorf("Error in TestVcfAnnotateAncestorFromFa. Expected: %s. Found: %s.", dna.BasesToString(answers[i]), dna.BasesToString(VcfQueryAncestor(v)))
		}
		i++
	}
}
