package bed

import (
	"log"
	"strings"
)

//Helper function to get regions as a bed that comes prev this bed
func negatePrevBed(curr *Bed, prev int) (int, *Bed) {
	if curr.ChromStart == 0 {
		log.Fatalf("Error: this helper function does not deal with the 0th base...\n")
		return 0, nil
	} else {
		return int(curr.ChromEnd), &Bed{Chrom: curr.Chrom, ChromStart: int64(prev), ChromEnd: curr.ChromStart}
	}
}

func InvertRegions(beds []*Bed, chromLen int) []*Bed {
	var ans []*Bed
	if len(beds) < 1 {
		log.Fatalf("Error: bed slice needs to contain at least one bed record...\n")
	} else if strings.Compare(beds[0].Chrom, beds[len(beds)-1].Chrom) != 0 {
		log.Fatalf("Error: chromosome of first bed in slice must be the same as the last slice...\n")
	} else {
		var prev int = 0
		var curr *Bed
		if beds[0].ChromStart == 0 {
			prev = int(beds[0].ChromEnd)
		} //else {
		//	ans = append(ans, &Bed{Chrom: beds[0].Chrom, ChromStart: 0, ChromEnd: beds[0].ChromEnd})

		//}
		for i := 0; i < len(beds); i++ {
			prev, curr = negatePrevBed(beds[i], prev)
			ans = append(ans, curr)
		}
		//Once we end loop, we deal with the right side
		if prev == chromLen {
			return ans
		} else {
			ans = append(ans, &Bed{Chrom: beds[len(beds)-1].Chrom, ChromStart: int64(prev), ChromEnd: int64(chromLen)})
		}
	}
	return ans
}
