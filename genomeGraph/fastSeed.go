package genomeGraph

import (
	"log"
)

func printSeedDev(a []*SeedDev) {
	for i := range a {
		log.Printf("tId:%d\ttStart:%d\tqStart:%d\tLen:%d\ttotalLen:%d\tposStrand:%t\n", a[i].TargetId, a[i].TargetStart, a[i].QueryStart, a[i].Length, a[i].TotalLength, a[i].PosStrand)
	}
}

func SortSeedDevByLen(seeds []*SeedDev) {
	stableSort(seeds, func(i, j int) bool { return CompareLenSeedDev(seeds[i], seeds[j]) == 1 })
}

func stableSort(a []*SeedDev, less func(i, j int) bool) {
	n := len(a)
	aux := make([]*SeedDev, n)
	copy(aux, a)
	mergeSort(a, aux, 0, n, less)
}

func mergeSort(a, aux []*SeedDev, lo, hi int, less func(i, j int) bool) {
	if hi-lo <= 1 {
		return
	}
	mid := lo + (hi-lo)/2
	mergeSort(aux, a, lo, mid, less)
	mergeSort(aux, a, mid, hi, less)
	merge(a, aux, lo, mid, hi, less)
}

func merge(a, aux []*SeedDev, lo, mid, hi int, less func(i, j int) bool) {
	i, j := lo, mid
	for k := lo; k < hi; k++ {
		if i == mid {
			aux[k] = a[j]
			j++
		} else if j == hi {
			aux[k] = a[i]
			i++
		} else if less(j, i) {
			aux[k] = a[j]
			j++
		} else {
			aux[k] = a[i]
			i++
		}
	}
}

func CompareLenSeedDev(a *SeedDev, b *SeedDev) int {
	if a.TotalLength == b.TotalLength {
		return 0
	} else if a.TotalLength < b.TotalLength {
		return -1
	} else if a.TotalLength > b.TotalLength {
		return 1
	} else {
		log.Fatalf("Error: SeedDev total length compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.TotalLength, b.TargetId, b.TargetStart, b.TotalLength)
		return 0
	}
}

func CompareSeedDev(a *SeedDev, b *SeedDev) int {
	if a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.Length == b.Length {
		return 0
	} else if a.TargetId < b.TargetId ||
		(a.TargetId == b.TargetId && a.TargetStart < b.TargetStart) ||
		(a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.Length < b.Length) {
		return -1
	} else if a.TargetId > b.TargetId ||
		(a.TargetId == b.TargetId && a.TargetStart > b.TargetStart) ||
		(a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.Length > b.Length) {
		return 1
	} else {
		log.Fatalf("Error: SeedDev compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}
