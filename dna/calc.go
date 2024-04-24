package dna

import (
	"fmt"
	"log"
	"math"
)

// 0 is deltaH, 1 is deltaS
var nnTable = map[string][]float64{
	"selfComp": {0, -1.4},
	"endAT":    {2.3, 4.1},
	"endGC":    {0.1, -2.8},
	"AA":       {-7.9, -22.2},
	"TT":       {-7.9, -22.2},
	"AT":       {-7.2, -20.4},
	"TA":       {-7.2, -21.2},
	"CA":       {-8.5, -22.7},
	"TG":       {-8.5, -22.7},
	"GT":       {-8.4, -22.4},
	"AC":       {-8.4, -22.4},
	"CT":       {-7.8, -21.0},
	"AG":       {-7.8, -21.0},
	"GA":       {-8.2, -22.2},
	"TC":       {-8.2, -22.2},
	"CG":       {-10.6, -27.2},
	"GC":       {-9.8, -24.4},
	"GG":       {-8.0, -19.9},
	"CC":       {-8.0, -19.9},
}

// R is the universal gas constant
var R float64 = 1.987

// k is related to the concentration
var k float64 = (25 - (25 / 2)) * 1e-9

func identity(b Base, slc []int) {
	switch b {
	case 0, 3:
		slc[0]++
	case 1, 2:
		slc[1]++
	default:
		log.Fatalf("non ACTG bases are not currently supported for Tm calculations")
	}
}

// {A/T, G/C}
func getEnds(seq []Base) []int {
	var ans []int = []int{0, 0}
	identity(seq[0], ans)
	identity(seq[len(seq)-1], ans)
	return ans
}

func MeltingTemp(seq []Base) float64 {
	var selfComp bool
	var deltaS, deltaH float64
	var val []float64

	AllToUpper(seq)

	if CompareSeqsIgnoreCaseAndGaps(seq, ReverseComplementAndCopy(seq)) == 0 {
		selfComp = true
	} else {
		selfComp = false
	}

	endVals := getEnds(seq)
	deltaH += nnTable["endAT"][0] * float64(endVals[0])
	deltaS += nnTable["endAT"][1] * float64(endVals[0])
	deltaH += nnTable["endGC"][0] * float64(endVals[1])
	deltaS += nnTable["endGC"][1] * float64(endVals[1])

	seqString := BasesToString(seq)

	for i := 0; i < len(seqString)-1; i++ {
		val = nnTable[seqString[i:i+2]]
		deltaH += val[0]
		deltaS += val[1]
	}

	if selfComp {
		k = 25 * 1e-9
		deltaH += nnTable["selfComp"][0]
		deltaS += nnTable["selfComp"][1]
	}

	deltaS += 0.368 * float64(len(seq)-1) * math.Log(25)

	fmt.Println("deltaH: ", deltaH)
	fmt.Println("deltaS: ", deltaS)
	fmt.Println(R * math.Log(k))
	return (1000*deltaH)/(deltaS+(R*math.Log(k))) - 273.15
}
