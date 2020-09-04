package vcf

import (
	"testing"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	//DEBUG: "fmt"
	"fmt"
)

var AncestorAlleles []dna.Base = []dna.Base{dna.A, dna.T, dna.G, dna.C, dna.A, dna.A, dna.C}

func TestGVcfQueryAncestor(t *testing.T) {
	reader, _ := GoReadToChan("testdata/AncestorTest.vcf")
	var input []dna.Base
	var i int = 0
	var v *GVcf
	for each := range reader {
		v = VcfToGvcf(each)
		input = GVcfQueryAncestor(v)
		if input[0] != AncestorAlleles[i] {
			t.Errorf("Error in TestGVcfQueryAncestor. Input: %s. Expected: %s.", dna.BaseToString(input[0]), dna.BaseToString(AncestorAlleles[i]))
		}
		i++
	}
}

func TestGVcfAppendAncestor(t *testing.T) {
	reader, _ := GoReadToChan("testdata/Ancestor_No_Annotation.vcf")
	var i int = 0
	var input []dna.Base
	var allele []dna.Base = make([]dna.Base, 1)
	var v *GVcf

	for each := range reader {
		v = VcfToGvcf(each)
		allele[0] = AncestorAlleles[i]
		GVcfAppendAncestor(v, allele)
		
		input = GVcfQueryAncestor(v)
		if input[0] != AncestorAlleles[i] {
			t.Errorf("Error in TestGVcfAppendAncestor. Input: %s. Expected: %s.", dna.BaseToString(input[0]), dna.BaseToString(AncestorAlleles[i]))
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
	if after - before != 1 {
		t.Errorf("Error in TestAncestorFlagToHeader.")
	}

	//now we test a file that has no info columns in the header
	_, header = GoReadToChan("testdata/Ancestor_No_Info.vcf")
	before = len(header.Text)
	AncestorFlagToHeader(header)
	after = len(header.Text)
	//DEBUG: PrintHeader(header)
	if after - before != 1 {
		t.Errorf("Error in TestAncestorFlagToHeader.")
	}
}

var answers [][]dna.Base = [][]dna.Base{dna.StringToBases("A"), dna.StringToBases("T"), dna.StringToBases("A"), dna.StringToBases("CCT"), dna.StringToBases("N")}

func TestGVcfAnnotateAncestorFromFa(t *testing.T) {
	reader, _ := GoReadToChan("testdata/TestGVcfAnnotateAncestorFromFa.vcf")
	records := fasta.Read("testdata/testAncestorSequence.fa")
	var v *GVcf
	var i int = 0

	for each := range reader {
		v = VcfToGvcf(each)
		GVcfAnnotateAncestorFromFa(v, records)
		fmt.Printf("Answer: %s. Expected:%s. \n", dna.BasesToString(GVcfQueryAncestor(v)), dna.BasesToString(answers[i]))
		if !dna.SeqEqual(GVcfQueryAncestor(v), answers[i]) {
			t.Errorf("Error in TestGVcfAnnotateAncestorFromFa. Expected: %s. Found: %s.", dna.BasesToString(answers[i]), dna.BasesToString(GVcfQueryAncestor(v)))
		}
		i++
	}
}




