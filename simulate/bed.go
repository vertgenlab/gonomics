package simulate

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
)

// CountWindows counts thte total viable windows of length regionLength in the sequence searchSpace
func CountWindows(searchSpace []bed.Bed, regionLength int) (int, int) {
	var length, totalWindows int
	for i := 0; i < len(searchSpace); i++ {
		length = searchSpace[i].ChromEnd - searchSpace[i].ChromStart

		if length >= regionLength {
			totalWindows = totalWindows + (length - regionLength + 1)
		}
	}

	return length, totalWindows
}

// GoSimulateBed takes a searchSpace (represented by a noGap.bed input file, here as a parsed struct) and generates a
// number of regions (regionCount) of a specified length (regionLength) and sends the simulated regions to an output chan.
func GoSimulateBed(searchSpace []bed.Bed, regionCount int, regionLength int) <-chan bed.Bed {
	var tmp, chromWindows int
	
	length, totalWindows := CountWindows(searchSpace, regionLength)

	// this function generates new bed regions and sends them to a channel.
	c := make(chan bed.Bed, 1000)
	go func() {
		for i := 0; i < regionCount; i++ {
			c <- GenerateBedRegion(searchSpace, totalWindows, regionLength)
		}
		close(c)
	}()
	return c
}

// GenerateBedRegion searches the regions of searchSpace (a noGap.bed input file, as a parsed struct) and randomly selects a continuous region of length regionLength 
func GenerateBedRegion(searchSpace []bed.Bed, totalWindows int, regionLength int) {
	tmp = numbers.RandIntInRange(0, totalWindows)

	// iterating through each ungapped window
	for j := 0; j < len(searchSpace); j++ {
		Length = searchSpace[j].ChromEnd - searchSpace[j].ChromStart

		// check that the window can fit a sequence of that length e.g. want 10 bp, but window is 5bp
		chromWindows = Length - regionLength + 1

		if chromWindows < 1 {
			continues
		}
		if tmp-chromWindows > 0 {
			tmp = tmp - chromWindows
		} else {
			return bed.Bed{Chrom: searchSpace[j].Chrom, ChromStart: searchSpace[j].ChromStart + tmp - 1, ChromEnd: searchSpace[j].ChromStart + tmp - 1 + regionLength, Name: searchSpace[j]
			}
		}
	}
}