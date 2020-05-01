package alleles

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"testing"
)
// TODO: build more robust test
// TODO: build test for graph genome
func TestFindNewVariation(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	experimental := "testdata/human_chrM.sam"
	normal := "testdata/chrM_head.sam"
	answer := FindNewVariation(ref, experimental, normal, 0, 1, 20, 100)

	for data := range answer {
		vcf.PrintSingleLine(data)
	}
}
