package main

import (
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/numbers"
)

func findDist(bedpeIn string) {
	var normal float64
	var allNormal []float64
	records := bedpe.Read(bedpeIn)
	bedpe.SortByCoord(records)
	x, m, s := bedpe.FindStats(records)

	for a := range x {
		normal = numbers.NormalDist(x[a], m[a], s[a])
		allNormal = append(allNormal, normal)
	}
}
