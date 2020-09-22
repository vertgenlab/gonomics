package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
)
//SimulateAFS returns an allele frequency spectrum AFS struct for n individuals with k segregating sites
//with a selection parameter alpha.
/*func SimulateAFS(alpha float64, n int, k int) AFS {
	var answer AFS
	answer.sites = make([]*SegSite, 0)
	for i := 0; i < k; i++ {
		answer.sites = append(answer.sites, &SegSite{StationaritySample(n, alpha), n})
	}
	return answer
}*/

//StationaritySample returns an allele frequency i out of n individuals sampled from a stationarity 
//distribution with selection parameter alpha.
func StationaritySampler(alpha float64, samples int, maxSampleDepth int, bins int, xLeft float64, xRight float64, randSeed bool) []float64 {
	f := AFSStationarityClosure(alpha)
	return numbers.FastRejectionSampler(xLeft, xRight, f, bins, maxSampleDepth, samples, randSeed)
}
