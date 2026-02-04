package simulate

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// CountWindows counts the total viable windows of length regionLength in the sequence searchSpace
func CountWindows(searchSpace []bed.Bed, regionLength int) int {
	var length, totalWindows int
	for i := range searchSpace {
		length = searchSpace[i].ChromEnd - searchSpace[i].ChromStart

		if length >= regionLength {
			totalWindows += length - regionLength + 1
		}
	}

	return totalWindows
}

// GenerateBedRegion searches the regions of searchSpace (a noGap.bed input file, as a parsed struct) and randomly selects a continuous region of length regionLength. RandPos can be an int in [0, totalWindows-1] inclusive where totalWindows is the total viable windows of size regionLength in searchSpace.
func GenerateBedRegion(searchSpace []bed.Bed, randPos int, regionLength int) (bed.Bed, bool) {
	var chromWindows int
	var length int
	// iterating through each ungapped window
	for j := range searchSpace {
		length = searchSpace[j].ChromEnd - searchSpace[j].ChromStart

		// number of windows of length regionLength in the specific bed region
		// check that the window can fit a sequence of that length e.g. want 10 bp, but window is 5bp
		chromWindows = length - regionLength + 1

		if chromWindows < 1 {
			continue
		}

		// Decrement randomly generated overall position (corresponds to start of generated region) until it fits within a region
		// must have randPos < chromWindows (randPos is 0-indexed, chromWindows is not), at most randPos + 1 = chromWindows
		// e.g. randPos=0, chromWindows must be at least 1, 0-1 = -1
		if randPos-chromWindows > -1 {
			randPos -= chromWindows
		} else {
			if searchSpace[j].Name == "" {
				return bed.Bed{
				Chrom: searchSpace[j].Chrom, 
				ChromStart: searchSpace[j].ChromStart + randPos, 
				ChromEnd: searchSpace[j].ChromStart + randPos + regionLength, 
				FieldsInitialized: 3}, true
			}
			return bed.Bed{
				Chrom:             searchSpace[j].Chrom,
				ChromStart:        searchSpace[j].ChromStart + randPos,
				ChromEnd:          searchSpace[j].ChromStart + randPos + regionLength,
				Name:              searchSpace[j].Name,
				FieldsInitialized: 4}, true
		}
	}

	log.Panic("Unable to generate region")
	return bed.Bed{}, false
}

// GoSimulateBed takes a searchSpace (represented by a noGap.bed input file, here as a parsed struct) and generates a
// number of regions (regionCount) of a specified length (regionLength) and sends the simulated regions to an output chan.
func GoSimulateBed(searchSpace []bed.Bed, regionCount int, regionLength int) <-chan bed.Bed {
	var randPos int
	c := make(chan bed.Bed, 1000)

	totalWindows := CountWindows(searchSpace, regionLength)

	// this function generates new bed regions and sends them to a channel.
	go func() {
		for i := 0; i < regionCount; i++ {
			randPos = numbers.RandIntInRange(0, totalWindows)
			if region, found := GenerateBedRegion(searchSpace, randPos, regionLength); found {
				c <- region
			}
		}
		close(c)
	}()
	return c
}
