package qDna

import ()

type QBase struct {
	A float32
	C float32
	G float32
	T float32
}

type QFrag struct {
	Seq  []*QBase
	From []*Location
	Fwd  []*QAdj
	Rev  []*QAdj
}

type QAdj struct {
	Prob float32
	Next *QFrag
}

type Location struct {
	Assembly string
	Chr      string
	Start    int64
	End      int64
}
