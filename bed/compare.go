package bed

import (
	"errors"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"sort"
	"strings"
)

//SortByCoord sorts in place a slice of Bed structs by their genomic position.
func SortByCoord(bedFile []Bed) {
	sort.Slice(bedFile, func(i, j int) bool { return Compare(bedFile[i], bedFile[j]) == -1 })
}

//SortBySize sorts in place a slice of Bed structs by their size from low to high.
func SortBySize(bedFile []Bed) {
	sort.Slice(bedFile, func(i, j int) bool { return CompareSize(bedFile[i], bedFile[j]) == -1 })
}

//SortByChromEndByChrom sorts in place a slice of Bed structs by their chrom, and then by chromEnd, but not by chromStart.
func SortByChromEndByChrom(bedFile []Bed) {
	sort.Slice(bedFile, func(i, j int) bool { return CompareChromEndByChrom(bedFile[i], bedFile[j]) == -1 })
}

func MergeBeds(bedFile []Bed) []Bed {
	SortByCoord(bedFile)
	var i, j int
	for i = 0; i < len(bedFile)-1; {
		if !Overlap(bedFile[i], bedFile[i+1]) {
			i++
		} else {
			bedFile[i].ChromStart, bedFile[i].ChromEnd, bedFile[i].Score = numbers.Min(bedFile[i].ChromStart, bedFile[i+1].ChromStart), numbers.Max(bedFile[i].ChromEnd, bedFile[i+1].ChromEnd), bedFile[i].Score+bedFile[i+1].Score
			for j = i + 1; j < len(bedFile)-1; j++ {
				bedFile[j] = bedFile[j+1]
			}
			bedFile = bedFile[:len(bedFile)-1]
		}
	}
	return bedFile
}

//Adjacent returns true if two input Bed entries are adjacent (one immediately follows the other).
func Adjacent(alpha Bed, beta Bed) bool {
	if alpha.Chrom != beta.Chrom {
		return false
	}
	if alpha.ChromEnd == beta.ChromStart || alpha.ChromStart == beta.ChromEnd { //adjacency for these quantities is equality as bed coordinates are half-open (contains start, does not contain end).
		return true
	}
	return false
}

//Overlap returns true if two input Bed entries have an overlap of any kind.
func Overlap(alpha Bed, beta Bed) bool {
	if (numbers.Max(alpha.ChromStart, beta.ChromStart) < numbers.Min(alpha.ChromEnd, beta.ChromEnd)) && strings.Compare(alpha.Chrom, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

//OverlapCount returns the number of elements from list one that have any overlap with list two. Answers range from 0 to len(a).
//Input bed slices must be presorted with SortByCoord.
func OverlapCount(a []Bed, b []Bed) int {
	var count int = 0
	var aIndex, bIndex int

	for aIndex < len(a) && bIndex < len(b) {
		if Overlap(a[aIndex], b[bIndex]) {
			count++
			aIndex++
		} else if CompareChromEndByChrom(a[aIndex], b[bIndex]) < 0 {
			aIndex++
		} else {
			bIndex++
		}
	}
	return count
}

//OverlapLengthSum calculates the total number of overlapping bases between two sets of bed elements.
//Input bed slices must be presorted with SortByCoord
func OverlapLengthSum(a []Bed, b []Bed) int {
	var sum int = 0
	var aIndex, bIndex, oLen int
	for aIndex < len(a) && bIndex < len(b) {
		oLen = OverlapLength(a[aIndex], b[bIndex])
		if oLen != 0 {
			sum += oLen
		}
		if CompareChromEndByChrom(a[aIndex], b[bIndex]) < 0 {
			aIndex++
		} else {
			bIndex++
		}
	}
	return sum
}

//OverlapLength returns the number of bases for which two Bed entries overlap.
func OverlapLength(a Bed, b Bed) int {
	if !Overlap(a, b) {
		return 0
	}
	end := numbers.Min(a.ChromEnd, b.ChromEnd)
	start := numbers.Max(a.ChromStart, b.ChromStart)
	return end - start
}

//Compare returns zero for equal beds and otherwise returns the ordering of the two Bed entries. Used for SortByCoord.
func Compare(a Bed, b Bed) int {
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

//CompareSize returns zero for beds with an equal length (ChromEnd - ChromStart) and otherwise returns the ordering of the two Bed structs.
func CompareSize(a Bed, b Bed) int {
	sizeA := a.ChromEnd - a.ChromStart
	sizeB := b.ChromEnd - b.ChromStart
	if sizeA < sizeB {
		return -1
	}
	if sizeA > sizeB {
		return 1
	}
	return 0
}

//CompareChromEnd returns zero for beds with an equal ChromEnd position and otherwise returns the ordering of the two Bed entries by chromEnd.
func CompareChromEnd(a Bed, b Bed) int {
	if a.ChromEnd < b.ChromEnd {
		return -1
	}
	if a.ChromEnd > b.ChromEnd {
		return 1
	}
	return 0
}

//CompareChromEndByChrom compares beds by chromosome and then by chromEnd, but not by chromStart.
func CompareChromEndByChrom(a Bed, b Bed) int {
	chromComp := strings.Compare(a.Chrom, b.Chrom)
	if chromComp != 0 {
		return chromComp
	}
	return CompareChromEnd(a, b)
}

//CompareDistance compares beds by chromStart and chromEnd and returns the minimum distance between
//the two beds. return options; -1 (different chromosomes, no distance calculated),
//0 (overlap, minimum distance is 0), >=0 is the minimum distance.
func CompareDistance(a Bed, b Bed) int {
	if a.Chrom != b.Chrom {
		myError := errors.New("Could not calculate distance, different chromosomes")
		return myError
	}
	if Overlap(a, b) {
		return 0
	}
	if a.ChromStart-b.ChromEnd >= 0 { //only positive if bed "a" is downstream of "b" bed.
		return a.ChromStart - b.ChromEnd + 1
	}
	if b.ChromStart-a.ChromEnd >= 0 { //only positive if bed "b" is downstream of "a" bed.
		return b.ChromStart - a.ChromEnd + 1

	}
	log.Panic("something went wrong with CompareDistance")
	return -1 // maybe just return -1 here to indicate that we can't count a distance?
}

//AllAreEqual returns true if two input slices of Beds contain Bed entries that all return true for Equal.
func AllAreEqual(a []Bed, b []Bed) bool {
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

//Equal returns true if two input Bed entries have the same Chrom, ChromStart, and ChromEnd. False otherwise.
func Equal(a Bed, b Bed) bool {
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
