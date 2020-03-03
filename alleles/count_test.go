package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

func TestSamToAlleles(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	answer := SamToAlleles("testdata/human_chrM.sam", ref, 0)

	for value := range answer {
		fmt.Sprintln(value)
	}
}
