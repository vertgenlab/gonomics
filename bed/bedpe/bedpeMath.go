package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"log"
	"math"
)

//TODO:maybe to do the whole chrom at once I find the mode of the start positions and then calculate just for that spot?

// FindStats will return an ordered slice for each of amplitude, mean and standard deviation, where each position
// corresponds to the bin in the same order on the genome given from a straw file (see hic package) that has been turned into a bedpe
func FindStats(in []BedPe) (x []float64, mu []float64, sigma []float64) {
	var amp, mean, sd []float64
	var s []bed.Bed
	var ok bool
	var i, n, k int
	var keys []int
	var dist, c, ampRelativePos float64
	matrix := make(map[int][]bed.Bed, 0) //chrom start of bin1 to []bed corresponds to all data in bedpe.B

	for i = range in {
		s, ok = matrix[in[i].A.ChromStart]
		if !ok {
			matrix[in[i].A.ChromStart] = append(matrix[in[i].A.ChromStart], in[i].B)
			keys = append(keys, in[i].A.ChromStart)
		} else {
			for b := range s {
				if s[b].ChromStart == in[i].B.ChromStart {
					s[b].Score += in[i].B.Score
				} else {
					matrix[in[i].A.ChromStart] = append(matrix[in[i].A.ChromStart], in[i].B)
				}
			}
		}
	}

	for _, k = range keys {
		c = 0
		dist = 0
		var dists []float64
		for n = range matrix[k] {
			dist = float64(matrix[k][n].ChromStart-k) * float64(matrix[k][n].Score)
			dists = append(dists, dist)
			c += float64(matrix[k][n].Score) //number of contacts for denominator of mean calculation
		}
		mean = append(mean, dist/c)
		sd = append(sd, calculateSd(mean[len(mean)-1], dists))
		ampRelativePos = float64(matrix[k][n].ChromStart + k)
		amp = append(amp, findAmplitude(matrix, ampRelativePos, k))
	}

	if len(amp) != len(mean) || len(amp) != len(sd) || len(mean) != len(sd) || len(mean) != len(keys) {
		log.Fatalf("Amplitude, mean or standard deviation are not of the correct length.\n amp=%v \n mu=%v \n sd =%v\n All Should be length=%v", len(amp), len(mean), len(sd), len(keys))
	}

	return amp, mean, sd
}

// calculateSd will determine the standard deviation for a set of values
func calculateSd(mean float64, values []float64) float64 {
	var sd float64
	for s := range values {
		sd += math.Pow(values[s]-mean, 2)
	}
	sd = math.Sqrt(sd / float64(len(values)))
	return sd
}

// findAmplitude will determine the score of the bedpe which contains the mean distance from thisBin
func findAmplitude(m map[int][]bed.Bed, relativePos float64, thisBin int) float64 {
	var answer float64

	for i := range m[thisBin] {
		if i == len(m[thisBin])-1 {
			return float64(m[thisBin][i].Score)
		} else {
			if thisBin-m[thisBin][i].ChromStart <= int(relativePos) && int(relativePos) <= thisBin-m[thisBin][i+1].ChromStart {
				return float64(m[thisBin][i].Score)
			}
		}
	}

	return answer
}
