package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

//TODO: build more robust test
func TestSyncAlleleStreams(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	one := SamToAlleles("testdata/chrM_head.sam", ref, 0)
	two := SamToAlleles("testdata/chrM_tail.sam", ref, 0)

	answer := SyncAlleleStreams(ref, 1000, one, two)

	for i := range answer {
		fmt.Sprintln("Finished", i[0].Location, len(i))
	}
}
