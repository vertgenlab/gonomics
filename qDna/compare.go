package qDna

/*
import (
	"log"
	"sort"
	"strings"
)

func CompareCoord(alpha *NoGapAln, beta *NoGapAln) int {
	if alpha.Start < beta.Start {
		return -1
	}
	if alpha.Start > beta.Start {
		return 1
	}
	return 0
}

func Sort(seedList []*NoGapAln) {
	sort.Slice(seedList, func(i, j int) bool { return CompareCoord(seedList[i], seedList[j]) == -1 })
}

func mergeTwoSeed(alpha *NoGapAln, beta *NoGapAln) []*NoGapAln {
	var answer []*NoGapAln
	var curr *NoGapAln
	if beta.Start <= alpha.Start && beta.End <= alpha.End && beta.End >= alpha.Start {
		curr = &NoGapAln{Start: beta.Start, End: alpha.End, Score: 0}
		answer = append(answer, curr)
	} else if beta.Start >= alpha.Start && beta.Start <= alpha.End && alpha.End <= beta.End {
		curr = &NoGapAln{Start: alpha.Start, End: beta.End, Score: 0}
		answer = append(answer, curr)
	} else if beta.Start >= alpha.Start && beta.Start <= alpha.End && beta.End >= alpha.Start && beta.End <= alpha.End {
		curr = &NoGapAln{Start: alpha.Start, End: alpha.End, Score: 0}
		answer = append(answer, curr)
	} else if alpha.Start >= beta.Start && alpha.Start <= beta.End && alpha.End >= beta.Start && alpha.End <= beta.End {
		curr = &NoGapAln{Start: beta.Start, End: beta.End, Score: 0}
		answer = append(answer, curr)
	} else {
		answer = append(answer, alpha)
		answer = append(answer, beta)
	}
	return answer
}

func MergedAllSeedRegions(seedList []*NoGapAln) []*NoGapAln {
	var answer []*NoGapAln
	Sort(allSeedRegions)
	if len(seedList) == 0 {
			log.Fatalf("Seed regions list is empty")
	}
	alpha := seedList[0]
	var beta *NoGapAln
	var curr *NoGapAln
	for i := 1; i < len(seedList); i++ {
		beta = seedList[i]
		if beta.Start <= alpha.End && beta.End <= alpha.End {
			curr = &NoGapAln{Start: alpha.Start, End: alpha.End, Score: 0}
		} else if beta.Start <= alpha.End && beta.End >= alpha.End {
			curr = &NoGapAln{Start: alpha.Start, End: beta.End, Score: 0}
		} else if beta.Start > alpha.End {

		} else {
			answer = append(answer, alpha)
			alpha = beta
		}







	}*/
