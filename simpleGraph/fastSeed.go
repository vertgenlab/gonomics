package simpleGraph

import (
	"log"
	"sort"
)

func printSeedDev(a []*SeedDev) {
	for i := range a {
		log.Printf("tId:%d\ttStart:%d\tqStart:%d\tLen:%d\ttotalLen:%d\tposStrand:%t\n", a[i].TargetId, a[i].TargetStart, a[i].QueryStart, a[i].Length, a[i].TotalLength, a[i].PosStrand)
	}
}

func SortSeedDevByLen(seeds []*SeedDev) {
	sort.Slice(seeds, func(i, j int) bool { return CompareLenSeedDev(seeds[i], seeds[j]) == 1 })
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
