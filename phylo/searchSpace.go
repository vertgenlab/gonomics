package phylo

import (
	"github.com/vertgenlab/gonomics/bed"
)

// MakeBitArrayFromSearchSpaceBed.
func MakeBitArrayFromSearchSpaceBed(searchSpaceFile string, referenceLength int, Chrom string) []bool {
	searchSpace := bed.Read(searchSpaceFile)
	bitArray := make([]bool, referenceLength)
	var j int

	for i := range searchSpace {
		if searchSpace[i].Chrom == Chrom {
			for j = searchSpace[i].ChromStart; j < searchSpace[i].ChromEnd; j++ {
				bitArray[j] = true
			}
		}
	}
	return bitArray
}
