package genomeGraph

import (
	"log"
	"sort"

	"github.com/vertgenlab/gonomics/fastq"
)

func getLastPart(a *SeedDev) *SeedDev {
	for a.NextPart != nil {
		a = a.NextPart
	}
	return a
}

func CompareBlastScore(a *SeedDev, b *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int {
	if BlastSeed(a, read, scoreMatrix) == BlastSeed(b, read, scoreMatrix) {
		return 0
	} else if BlastSeed(a, read, scoreMatrix) < BlastSeed(b, read, scoreMatrix) {
		return -1
	} else if BlastSeed(a, read, scoreMatrix) > BlastSeed(b, read, scoreMatrix) {
		return 1
	} else {
		log.Fatalf("Error: SeedDev len compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}

func SortSeedsByScore(seeds []*SeedDev, read fastq.Fastq, scoreMatrix [][]int64) {
	sort.Slice(seeds, func(i, j int) bool {
		return BlastSeed(seeds[i], read, scoreMatrix) > BlastSeed(seeds[j], read, scoreMatrix)
	})
}

func SortBlastz(seeds []*SeedDev, read fastq.Fastq, scoreMatrix [][]int64) {
	sort.Slice(seeds, func(i, j int) bool { return CompareBlastScore(seeds[i], seeds[j], read, scoreMatrix) == 1 })
}

func BlastSeed(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	if seed.NextPart == nil {
		return scoreSeed(seed, read, scoreMatrix)
	} else {
		return scoreSeed(seed, read, scoreMatrix) + scoreSeed(seed.NextPart, read, scoreMatrix)
	}
}
