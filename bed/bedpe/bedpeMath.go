package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"math"
	"sort"
)

func FindStats(in []BedPe) (x []float64, mu []float64, sigma []float64) {
	var amp, mean, sd []float64
	var s []bed.Bed
	var ok bool
	var i, n, k int
	var keys []int
	var dist, c, ampRelativePos float64
	var matrix map[int][]bed.Bed //chrom start of bin to []bed
	//TODO: sort bedpe by start of bed A

	for i = range in {
		s, ok = matrix[in[i].A.ChromStart]
		if !ok {
			s = append(s, in[i].B)
			keys = append(keys, in[i].A.ChromStart)
		} else {
			for b := range s {
				if s[b].ChromStart == in[i].B.ChromStart {
					s[b].Score += in[i].B.Score
				}
			}
		}
	}

	sort.Ints(keys)
	//find mean sd and amplitude for each bin contacted by each bin
	for k = range keys {
		c = 0
		dist = 0
		var dists []float64
		for n = range matrix[k] {
			//score += float64(matrix[k][n].Score)
			dist = float64(matrix[k][n].ChromStart-k) * float64(matrix[k][n].Score)
			dists = append(dists, dist)
			c += float64(matrix[k][n].Score) //number of contacts for denominator to mean
		}
		mean = append(mean, dist/c)
		sd = append(sd, calculateSd(mean[len(mean)-1], dists))
		ampRelativePos = float64(matrix[k][n].ChromStart + k)
		amp = append(amp, findAmplitude(matrix, keys, ampRelativePos, k))
	}

	return amp, mean, sd
}

func calculateSd(mean float64, values []float64) float64 {
	var sd float64
	for s := range values {
		sd += math.Pow(values[s]-mean, 2)
	}
	sd = math.Sqrt(sd / float64(len(values)))
	return sd
}

func findAmplitude(m map[int][]bed.Bed, orderedKeys []int, relativePos float64, thisBin int) float64 {
	var answer float64

	for i := range m[thisBin] {
		if thisBin-m[thisBin][i].ChromStart <= int(relativePos) && int(relativePos) <= thisBin-m[thisBin][i+1].ChromStart {
			return float64(m[thisBin][i].Score)
		}
	}

	return answer
}
