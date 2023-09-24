package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/interval"
)

func main() {
	small := bed.Read("testdata/smallAJ.bed")
	large := bed.Read("testdata/largeAJ.bed")
	answer := asymetricJaccard(small, large, 1)
	fmt.Println("asymetric Jaccard: ", answer)
	answer = kulczynski2(small, large)
	fmt.Println("kulczynski2: ", answer)
	answer = kulczynski2point0(small, large)
	fmt.Println("kulczynski2.0: ", answer)
}

func asymetricJaccard(small []bed.Bed, big []bed.Bed, smallWindowSize int) float64 {
	var overlap []interval.Interval
	var totalOverlap int
	var nonOverlaps int
	var currOverlap interval.Interval

	smallIntervalMap := makeIntervalMap(small)

	for _, bd := range big {
		overlap = interval.Query(smallIntervalMap, bd, "any")
		if len(overlap) == 0 {
			nonOverlaps++
		}
		for _, currOverlap = range overlap {
			totalOverlap += overlapSize(bd, currOverlap.(bed.Bed))
		}
	}
	totalSizeSmall := bed.TotalSize(small)
	return float64(totalOverlap) / (float64(nonOverlaps*smallWindowSize) + float64(totalSizeSmall))
}

func kulczynski2(small []bed.Bed, large []bed.Bed) float64 {
	// 0.5 * (a/(a+b) + a/(a+c))
	// a = present in both samples
	// b = only in small sample
	// c = only in large sample

	var overlap []interval.Interval
	var a, b, c int

	smallIntervalMap := makeIntervalMap(small)

	for _, bd := range large {
		overlap = interval.Query(smallIntervalMap, bd, "any")
		a += len(overlap)
		if len(overlap) == 0 {
			c++
		} else {
			//a++
		}
	}

	bigIntervalMap := makeIntervalMap(large)
	for _, bd := range small {
		overlap = interval.Query(bigIntervalMap, bd, "any")
		if len(overlap) == 0 {
			b++
		}
	}
	return 0.5 * ((float64(a) / float64(a+b)) + (float64(a) / float64(a+c)))
}

func kulczynski2point0(small []bed.Bed, large []bed.Bed) float64 {
	var overlapSmall, overlapLarge, allLargeOverlaps, allSmallOverlaps []interval.Interval

	smallIntervalMap := makeIntervalMap(small)
	largeIntervalMap := makeIntervalMap(large)
	for _, bd := range large {
		overlapSmall = interval.Query(smallIntervalMap, bd, "any")
		for _, i := range overlapSmall {
			allSmallOverlaps = append(allSmallOverlaps, i)
		}
	}
	for _, bd := range small {
		overlapLarge = interval.Query(largeIntervalMap, bd, "any")
		for _, i := range overlapLarge {
			allLargeOverlaps = append(allLargeOverlaps, i)
		}
	}

	allSmallOverlapsBed := sliceIntervalToBed(allSmallOverlaps)
	allLargeOverlapsBed := sliceIntervalToBed(allLargeOverlaps)
	bed.SortByCoord(allSmallOverlapsBed)
	bed.SortByCoord(allLargeOverlapsBed)
	uniqueSmallOverlapBeds := unique(allSmallOverlapsBed)
	uniqueLargeOverlapBeds := unique(allLargeOverlapsBed)

	smallOverlapPerc := float64(len(uniqueSmallOverlapBeds)) / float64(len(small))
	bigOverlapPerc := float64(len(uniqueLargeOverlapBeds)) / float64(len(large))

	return (smallOverlapPerc + bigOverlapPerc) / 2
}

func makeIntervalMap(inBed []bed.Bed) map[string]*interval.IntervalNode {
	var smallWindowIntervals []interval.Interval

	for bd := range inBed {
		smallWindowIntervals = append(smallWindowIntervals, inBed[bd])
	}
	return interval.BuildTree(smallWindowIntervals)
}

func overlapSize(a bed.Bed, b bed.Bed) int {
	if !bed.Overlap(a, b) {
		return 0
	}
	switch {
	case bed.Equal(a, b):
		return a.ChromEnd - a.ChromStart
	case a.ChromStart <= b.ChromStart && a.ChromEnd >= b.ChromEnd:
		return b.ChromEnd - b.ChromStart
	case a.ChromStart >= b.ChromStart && a.ChromEnd <= b.ChromEnd:
		return a.ChromEnd - a.ChromStart
	case a.ChromStart >= b.ChromStart && a.ChromEnd >= b.ChromEnd:
		return b.ChromEnd - a.ChromStart
	case a.ChromStart <= b.ChromStart && a.ChromEnd <= b.ChromEnd:
		return a.ChromEnd - b.ChromStart
	}
	return 0
}

func sliceIntervalToBed(intervals []interval.Interval) []bed.Bed {
	var bedSlice []bed.Bed
	for _, i := range intervals {
		bedSlice = append(bedSlice, i.(bed.Bed))
	}
	return bedSlice
}

func unique(regions []bed.Bed) []bed.Bed {
	var uniqueRegions []bed.Bed

	currRegion := regions[0]
	for _, i := range regions {
		if bed.Equal(currRegion, i) {
			continue
		}
		uniqueRegions = append(uniqueRegions, currRegion)
		currRegion = i
	}
	uniqueRegions = append(uniqueRegions, currRegion)
	return uniqueRegions
}
