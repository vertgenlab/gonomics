package main

import "github.com/vertgenlab/gonomics/sam"

// CacheStruct contains a collection of all caches used for the likelihood and
// prior calculations in samAssembler.
type CacheStruct struct {
	DiploidBasePriorCache  [][]float64
	DiploidIndelPriorCache []float64
	HaploidBasePriorCache  [][]float64
	HaploidIndelPriorCache []float64
	HomozygousBaseCache    [][]float64
	HeterozygousBaseCache  [][]float64
	HomozygousIndelCache   [][]float64
	HeterozygousIndelCache [][]float64
}

// cacheSetup returns a CacheStruct of initialized samAssembler caches based on the user Settings struct.
func cacheSetup(s Settings) CacheStruct {
	var i int
	var diploidBasePriorCache [][]float64
	if s.FlatPrior {
		diploidBasePriorCache = sam.MakeDiploidBaseFlatPriorCache()
	} else {
		diploidBasePriorCache = sam.MakeDiploidBasePriorCache(s.Delta, s.Gamma)
	}
	var diploidIndelPriorCache []float64 = sam.MakeDiploidIndelPriorCache(s.Kappa, s.Delta)
	var haploidBasePriorCache [][]float64 = sam.MakeHaploidBasePriorCache(s.Delta, s.Gamma)
	var haploidIndelPriorCache []float64 = sam.MakeHaploidIndelPriorCache(s.Delta, s.Kappa)
	var homozygousBaseCache [][]float64 = make([][]float64, s.LikelihoodCacheSize)
	for i = range homozygousBaseCache {
		homozygousBaseCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	var heterozygousBaseCache = make([][]float64, s.LikelihoodCacheSize)
	for i = range heterozygousBaseCache {
		heterozygousBaseCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	var homozygousIndelCache = make([][]float64, s.LikelihoodCacheSize)
	for i = range homozygousIndelCache {
		homozygousIndelCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	var heterozygousIndelCache = make([][]float64, s.LikelihoodCacheSize)
	for i = range heterozygousIndelCache {
		heterozygousIndelCache[i] = make([]float64, s.LikelihoodCacheSize)
	}

	return CacheStruct{
		DiploidBasePriorCache:  diploidBasePriorCache,
		DiploidIndelPriorCache: diploidIndelPriorCache,
		HaploidBasePriorCache:  haploidBasePriorCache,
		HaploidIndelPriorCache: haploidIndelPriorCache,
		HomozygousBaseCache:    homozygousBaseCache,
		HeterozygousBaseCache:  heterozygousBaseCache,
		HomozygousIndelCache:   homozygousIndelCache,
		HeterozygousIndelCache: heterozygousIndelCache,
	}
}
