package simulate

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
)

// GoSimulateBed takes a searchSpace (represented by a noGap.bed input file, here as a parsed struct) and generates a
// number of regions (regionCount) of a specified length (regionLength) and sends the simulated regions to an output chan.
func GoSimulateBed(searchSpace []bed.Bed, regionCount int, regionLength int) <-chan bed.Bed {
	var Length, tmp, chromWindows int
	var totalWindows int
	c := make(chan bed.Bed, 1000)

	//count total viable windows
	for i := 0; i < len(searchSpace); i++ {
		Length = searchSpace[i].ChromEnd - searchSpace[i].ChromStart

		if Length >= regionLength {
			totalWindows = totalWindows + (Length - regionLength + 1)
		}
	}

	//this function generates new bed regions and sends them to a channel.
	go func() {
		for i := 0; i < regionCount; i++ {
			tmp = numbers.RandIntInRange(0, totalWindows)
			for j := 0; j < len(searchSpace); j++ {
				Length = searchSpace[j].ChromEnd - searchSpace[j].ChromStart
				chromWindows = Length - regionLength + 1
				//is chrom big enough?
				if chromWindows < 1 {
					continue
				}
				if tmp-chromWindows > 0 {
					tmp = tmp - chromWindows
				} else {
					c <- bed.Bed{Chrom: searchSpace[j].Chrom, ChromStart: searchSpace[j].ChromStart + tmp - 1, ChromEnd: searchSpace[j].ChromStart + tmp - 1 + regionLength, Name: searchSpace[j].Name, FieldsInitialized: 4}
					break
				}
			}
		}
		close(c)
	}()
	return c
}
