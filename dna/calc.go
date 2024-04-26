package dna

import (
	"log"
	"math"
)

// nnTable is a map that has nearest neighbor dinucleotides and they're corresponding deltaH () and deltaS () values, respectively.
// Entropy/enthalpy values were empirically determined in: Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594
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

// R is the universal gas constant (Cal/degrees C*Mol)
var R float64 = 1.987

// k is related to the concentration. Assumes both the piece of DNA and what it's annealing to are both at 50 nM
var k float64 = (25 - (25 / 2)) * 1e-9

// identity is a helper function for MeltingTemp. It takes in a base and slice of int with length 2. The function iterates
// the value in slice[0] if the base is an A/T or iterates slice[1] if the base is a C/G. Log fatal if non-ACTG bases are provided
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

// getEnds is a helper function for MeltingTemp that returns the number of AT or CG bases at the ends of a slice of Base.
// The returned slice has number of AT bases at index 0, and number of GC bases at index 1. Log fatals if non-ACTG bases are present at the terminal positions
func getEnds(seq []Base) []int {
	var ans []int = []int{0, 0}
	identity(seq[0], ans)
	identity(seq[len(seq)-1], ans)
	return ans
}

// MeltingTemp calculates the melting temp of slice of Base in Celsius with the nearest-neighbor algorithm. Assumes 50nM of
// both oligo and template and 50 mM Na+.
func MeltingTemp(seq []Base) float64 {
	var selfComp bool
	var deltaS, deltaH float64
	var val []float64

	AllToUpper(seq)

	//check to see if the oligo sequence is perfectly self-complimentary
	if CompareSeqsIgnoreCaseAndGaps(seq, ReverseComplementAndCopy(seq)) == 0 {
		selfComp = true
	} else {
		selfComp = false
	}

	//get starting deltaS and deltaH values based on nucleation bases
	endVals := getEnds(seq)
	deltaH += nnTable["endAT"][0] * float64(endVals[0])
	deltaS += nnTable["endAT"][1] * float64(endVals[0])
	deltaH += nnTable["endGC"][0] * float64(endVals[1])
	deltaS += nnTable["endGC"][1] * float64(endVals[1])

	seqString := BasesToString(seq)
	//iterate over all dinucleotides in the provided sequence
	for i := 0; i < len(seqString)-1; i++ {
		val = nnTable[seqString[i:i+2]]
		deltaH += val[0]
		deltaS += val[1]
	}

	//if the sequences are fully self-complimentary, a special correction needs to be applied
	if selfComp {
		k = 25 * 1e-9
		deltaH += nnTable["selfComp"][0]
		deltaS += nnTable["selfComp"][1]
	}

	//salt concentration correction
	deltaS += 0.368 * float64(len(seq)-1) * math.Log(50*1e-3)

	//full melting temp formula with all the parts. (-273.15 is the K to C conversion)
	return (1000*deltaH)/(deltaS+(R*math.Log(k))) - 273.15
}
