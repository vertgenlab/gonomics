package sam

import (
	"testing"
)

func TestSeekBamRegion(t *testing.T) {
	br, _ := OpenBam("testdata/rand.bam")
	bai := ReadBai("testdata/rand.bam.bai")
	reads := SeekBamRegion(br, bai, "chr7", 45000000, 45200000)
	for i := range reads {
		if reads[i].RName != "chr7" || !(reads[i].GetChromStart() < 45200000 && reads[i].GetChromEnd() > 45000000) {
			t.Error("problem with SeekBamRegion")
		}
	}
}
