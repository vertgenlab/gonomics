package simulate

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
	"sync"
)

func GoSimulateBed(searchSpace []*bed.Bed, regionCount int, regionLength int64) <-chan bed.Bed {
	var Length, tmp, chromWindows int64
	var totalWindows int
	c := make(chan bed.Bed)

	//count total viable windows
	for i := 0; i < len(searchSpace); i++ {
		Length = searchSpace[i].ChromEnd - searchSpace[i].ChromStart

		if Length >= regionLength {
			totalWindows = totalWindows + int(Length-regionLength)
		}
	}

	var wg sync.WaitGroup
	wg.Add(1)

	// DEBUG fmt.Printf("totalWindows: %d\n", totalWindows)
	go func() {
		for i := 0; i < regionCount; i++ {
				tmp = int64(numbers.RandIntInRange(0, totalWindows))
				// DEBUG fmt.Printf("Random number is %d\n", tmp)
				for j := 0; j < len(searchSpace); j++ {
					// DEBUG fmt.Printf("ChromEnd:%v ChromStart:%v \n",noGap[j].ChromEnd,noGap[j].ChromStart)
					Length = searchSpace[j].ChromEnd - searchSpace[j].ChromStart
					// DEBUG fmt.Printf("Length; %v.\n", Length)
					chromWindows = Length - regionLength + 1
					// DEBUG fmt.Printf("j; %v.\n", j)
					//is chrom big enough?
					if chromWindows < 1 {
						continue
					}
					if tmp-chromWindows > 0 {
						tmp = tmp - chromWindows
					} else {
						// DEBUG fmt.Printf("Got one\n")
						c <- bed.Bed{Chrom: searchSpace[j].Chrom, ChromStart: searchSpace[j].ChromStart + tmp - 1, ChromEnd: searchSpace[j].ChromStart + tmp - 1 + regionLength, Name: searchSpace[j].Name}
						break
					}
				}
			}
		wg.Done()
	}()
	go func() {
		wg.Wait()
		close(c)
	}()
	return c
}
