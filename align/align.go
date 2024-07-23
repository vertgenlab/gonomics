// Package align provides basic struct(s), gap scoring matrices, and dynamical approaches to aligning sequences
package align

import (
	"math"

	"github.com/vertgenlab/gonomics/cigar"
)

var veryNegNum int64 = math.MinInt64 / 2

// DefaultScoreMatrix is a DNA-DNA scoring matrix that pairs well with opening and extension penalites of:
// O=400 E=30.  It may be good for distances similar to human-mouse.
var DefaultScoreMatrix = [][]int64{
	{91, -114, -31, -123, -44},
	{-114, 100, -125, -31, -43},
	{-31, -125, 100, -114, -43},
	{-123, -31, -114, 91, -44},
	{-44, -43, -43, -44, -43},
}

// DefaultScoreMatrix is a DNA-DNA scoring matrix that pairs well with opening and extension penalites of:
// O=400 E=30.  It may be good for distances similar to human-fish.
var HoxD55ScoreMatrix = [][]int64{
	{91, -114, -31, -123, 0},
	{-114, 100, -125, -31, 0},
	{-31, -125, 100, -114, 0},
	{-123, -31, -114, 91, 0},
	{0, 0, 0, 0, 0},
}

// DefaultScoreMatrix is a DNA-DNA scoring matrix that pairs well with opening and extension penalites of:
// O=600 E=55.  It may be good for distances similar to mouse-rat.
var MouseRatScoreMatrix = [][]int64{
	{91, -114, -31, -123, 0},
	{-114, 100, -125, -31, 0},
	{-31, -125, 100, -114, 0},
	{-123, -31, -114, 91, 0},
	{0, 0, 0, 0, 0},
}

// DefaultScoreMatrix is a DNA-DNA scoring matrix that pairs well with opening and extension penalites of:
// O=600 E=150.  It may be good for distances similar to human-chimp.
var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

/*
var StrictScoreMatrix = [][]int64{
        {91, -114, -31, -123, -44},
        {-114, 100, -125, -31, -43},
        {-31, -125, 100, -114, -43},
        {-123, -31, -114, 91, -44},
        {-44, -43, -43, -44, -43},
}
*/


func countAlignmentColumns(route []cigar.Cigar) int64 {
	var count int = 0
	for i := range route {
		count += route[i].RunLength
	}
	return int64(count)
}
