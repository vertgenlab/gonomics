package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

func TestSyncAlleleStreams(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	one := SamToAlleles("testdata/human_chrM.sam", ref, 0)
	two := SamToAlleles("testdata/human_chrM.sam", ref, 0)

	answer := SyncAlleleStreams(ref, one, two)

	for i := range answer {
		fmt.Println("finished", i[0].Location, len(i))
	}
}