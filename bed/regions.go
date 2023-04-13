package bed

import (
	"log"
	"strings"
)

// negatePrevBed is a helper function to get regions as a bed that comes prev this bed
func negatePrevBed(curr Bed, prev int) (int, Bed) {
	if curr.ChromStart == 0 {
		log.Fatalf("Error: this helper function does not deal with the 0th base...\n")
		return 0, Bed{}
	} else {
		return int(curr.ChromEnd), Bed{Chrom: curr.Chrom, ChromStart: prev, ChromEnd: curr.ChromStart}
	}
}

// InvertRegions uses a set of bed regions, to return bed regions that are inverted, negated, and/or not covered by the given input. Bed entries should be sorted by chromosome position and cannot contain overlaps.
func InvertRegions(beds []Bed, chromLen int) []Bed {
	var ans []Bed
	if len(beds) < 1 {
		log.Fatalf("Error: bed slice needs to contain at least one bed record...\n")
	} else if strings.Compare(beds[0].Chrom, beds[len(beds)-1].Chrom) != 0 {
		log.Fatalf("Error: chromosome of first bed in slice must be the same as the last slice...\n")
	} else {
		var prev int = 0
		var i int = 0
		var curr Bed
		//TODO: Consider handling the start and end positions differently (not satisfied with this solution)
		if beds[0].ChromStart == 0 {
			i++
			prev = int(beds[0].ChromEnd)
		}
		for ; i < len(beds); i++ {
			prev, curr = negatePrevBed(beds[i], prev)
			ans = append(ans, curr)
		}
		//Once we end loop, we deal with the right side
		if prev == chromLen {
			return ans
		} else {
			ans = append(ans, Bed{Chrom: beds[len(beds)-1].Chrom, ChromStart: prev, ChromEnd: chromLen})
		}
	}
	return ans
}
