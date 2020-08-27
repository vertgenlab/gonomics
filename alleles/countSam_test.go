package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

func TestSamToAlleles(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	answer := GoCountSamAlleles("testdata/human_chrM.sam", ref, 0)

	var i int
	for value := range answer {
		if i < 10 {
			fmt.Sprintln(value.Count)
		}
		i++
	}
}

func TestNewCountAlleles(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	answerChan := GoCountSamAlleles("testdata/human_chrM.sam", ref, 0)

	var i int
	for value := range answerChan {
		if i < 10 {
			fmt.Sprintln(value.Count)
		}
		i++
	}
}
