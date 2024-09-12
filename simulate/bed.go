package simulate

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// CountWindows counts the total viable windows of length regionLength in the sequence searchSpace
func CountWindows(searchSpace []bed.Bed, regionLength int) int {
	var length, totalWindows int
	for i := 0; i < len(searchSpace); i++ {
		length = searchSpace[i].ChromEnd - searchSpace[i].ChromStart

		if length >= regionLength {
			totalWindows = totalWindows + (length - regionLength + 1)
		}
	}

	return totalWindows
}

// GenerateBedRegion searches the regions of searchSpace (a noGap.bed input file, as a parsed struct) and randomly selects a continuous region of length regionLength 
func GenerateBedRegion(searchSpace []bed.Bed, tmp int, regionLength int) (bed.Bed, bool) {
	var chromWindows int
	var length int
	// iterating through each ungapped window
	for j := 0; j < len(searchSpace); j++ {
		length = searchSpace[j].ChromEnd - searchSpace[j].ChromStart

		// check that the window can fit a sequence of that length e.g. want 10 bp, but window is 5bp
		chromWindows = length - regionLength + 1

		if chromWindows < 1 {
			continue
		}
		if tmp-chromWindows > 0 {
			log.Print("decrement\n")
			tmp -= chromWindows
		} else {
			log.Print("searchspace chromstart: ", searchSpace[j].ChromStart, "\t temp: ", tmp, "\n")
			log.Print("chromstart: ", searchSpace[j].ChromStart + tmp - 1, "\n")
			return bed.Bed{
				Chrom: searchSpace[j].Chrom, 
				ChromStart: searchSpace[j].ChromStart + tmp - 1, 
				ChromEnd: searchSpace[j].ChromStart + tmp - 1 + regionLength, 
				Name: searchSpace[j].Name,
				FieldsInitialized: 4}, true
		}
	}

	log.Panic("Error")
	return bed.Bed{}, false
}

// GoSimulateBed takes a searchSpace (represented by a noGap.bed input file, here as a parsed struct) and generates a
// number of regions (regionCount) of a specified length (regionLength) and sends the simulated regions to an output chan.
func GoSimulateBed(searchSpace []bed.Bed, regionCount int, regionLength int) <-chan bed.Bed {
	var tmp int
	c := make(chan bed.Bed, 1000)

	totalWindows := CountWindows(searchSpace, regionLength)

	//this function generates new bed regions and sends them to a channel.
	go func() {
		for i := 0; i < regionCount; i++ {
			tmp = numbers.RandIntInRange(0, totalWindows)
			if region, found := GenerateBedRegion(searchSpace, tmp, regionLength); found {
				c <- region
			}
		}
		close(c)
	}()
	return c
}
