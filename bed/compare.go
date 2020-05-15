package bed

import (
	"github.com/vertgenlab/gonomics/common"
	"sort"
	"strings"
	"log"
)

func SortByCoord(bedFile []*Bed) {
	sort.Slice(bedFile, func(i, j int) bool { return Compare(bedFile[i], bedFile[j]) == -1 })
}

func MergeBeds(bedFile []*Bed) []*Bed {
	SortByCoord(bedFile)
	var i, j int
	for i = 0; i < len(bedFile)-1; {
		if !Overlap(bedFile[i], bedFile[i+1]) {
			i++
		} else {
			bedFile[i].ChromStart, bedFile[i].ChromEnd, bedFile[i].Score = common.MinInt64(bedFile[i].ChromStart, bedFile[i+1].ChromStart), common.MaxInt64(bedFile[i].ChromEnd, bedFile[i+1].ChromEnd), bedFile[i].Score+bedFile[i+1].Score
			for j = i + 1; j < len(bedFile)-1; j++ {
				bedFile[j] = bedFile[j+1]
			}
			bedFile = bedFile[:len(bedFile)-1]
		}
	}
	return bedFile
}

func Overlap(alpha *Bed, beta *Bed) bool {
	if (common.MaxInt64(alpha.ChromStart, beta.ChromStart) < common.MinInt64(alpha.ChromEnd, beta.ChromEnd)) && strings.Compare(alpha.Chrom, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func OverlapLength(a *Bed, b *Bed) int64 {
	if !Overlap(a, b) {
		return 0
	}
	if a.ChromStart < b.ChromStart {
		return a.ChromEnd-b.ChromStart
	}else if a.ChromStart > b.ChromStart{
		return b.ChromEnd-a.ChromStart
	}
	if a.ChromEnd > b.ChromEnd {
		return b.ChromEnd-b.ChromStart
	}else if a.ChromEnd < b.ChromEnd {
			return a.ChromEnd-a.ChromStart
	}
	log.Fatalf("Critical Failure: End of Days")
	return -1
}

func Compare(a *Bed, b *Bed) int {
	chromComp := strings.Compare(a.Chrom, b.Chrom)
	if chromComp != 0 {
		return chromComp
	}
	if a.ChromStart < b.ChromStart {
		return -1
	}
	if a.ChromStart > b.ChromStart {
		return 1
	}
	if a.ChromEnd < b.ChromEnd {
		return -1
	}
	if a.ChromEnd > b.ChromEnd {
		return 1
	}
	return 0
}

func AllAreEqual(a []*Bed, b []*Bed) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !Equal(a[i], b[i]) {
			return false
		}
	}
	return true
}

func Equal(a *Bed, b *Bed) bool {
	if strings.Compare(a.Chrom, b.Chrom) != 0 {
		return false
	}
	if a.ChromStart != b.ChromStart {
		return false
	}
	if a.ChromEnd != b.ChromEnd {
		return false
	}
	return true
}
