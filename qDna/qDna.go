package qDna

import ()

type QBase struct {
	A float64
	C float64
	G float64
	T float64
}

type QFrag struct {
	Seq  []*QBase
	From []*Location
	Fwd  []*QAdj
	Rev  []*QAdj
}

type QAdj struct {
	Prob float64
	Next *QFrag
}

type Location struct {
	Assembly string
	Chr      string
	Start    int
	End      int
}
