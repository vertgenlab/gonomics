package maf

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta" 
	"testing"
	"log"
	"os"
)

var multiFaToMafTests = []struct {
	multifaFile 	string
	speciesGiven	bool
	speciesChr		[]string
	outDir			string
	expected		string
}{
	{"testdata/multifaToMaf_test_1.mfa", true, []string{"chr1", "chr2", "chr3", "chr4"}, "testdata/multifaToMaf_test_1.out.maf", "testdata/multifaToMaf_test_1.expected.maf"},
}

func TestMultiFaToMaf(t *testing.T) {
	var err error
	var out []*Maf
	var inputMultifa []fasta.Fasta
	for _, v := range multiFaToMafTests {
		inputMultifa = fasta.Read(v.multifaFile)
		out = MultifaToMaf(inputMultifa, v.speciesGiven, v.speciesChr)
		Write(v.outDir, out)
		log.Println(out)

		if !fileio.AreEqual(v.outDir, v.expected) {
			t.Errorf("Error in maf MultifaToMaf")
		} else {
			err = os.Remove(v.outDir)
		}
		exception.PanicOnErr(err)
	} 
}