package genomeGraph

import (
	"math"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

type Matrix struct {
	matrix [][]int64
	trace  [][]byte

	i       int
	j       int
	index   int
	currMax int64
	route   []cigar.Cigar
}

type Alignment struct {
	Seq         []dna.Base
	Path        []uint32
	queryStart  int
	queryEnd    int
	targetStart int
	targetEnd   int
	currScore   int64
}

type MemoryAllocation struct {
	Matrix    *Matrix
	Alignment *Alignment
}

type SeedMemory struct {
	currHits                                  []uint64
	codedNodeCoord                            uint64
	seqKey                                    uint64
	keyShift                                  uint
	keyIdx, keyOffset, readOffset, nodeOffset int
	nodeIdx, nodePos                          int64
	leftMatches, rightMatches                 int
	tempSeed                                  Seed
	Hits                                      []Seed
	Worker                                    []Seed
}

type scoreKeeper struct {
	targetStart  int
	targetEnd    int
	queryStart   int
	queryEnd     int
	queryLength  int
	currScore    int64
	seedScore    int64
	perfectScore int64
	leftScore    int64
	rightScore   int64
	leftPath     []uint32
	rightPath    []uint32
	leftSeq      []dna.Base
	rightSeq     []dna.Base
	currSeq      []dna.Base
	tailSeed     *Seed

	currSeed       *Seed
	leftAlignment  []cigar.Cigar
	rightAlignment []cigar.Cigar
}

func NewMemoryAllocation(size int) *sync.Pool {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return &sync.Pool{
		New: func() interface{} {
			pool := MemoryAllocation{
				Matrix: &Matrix{
					matrix:  m,
					trace:   trace,
					i:       0,
					j:       0,
					index:   0,
					currMax: math.MinInt64,
					route:   make([]cigar.Cigar, 0, 1),
				},
				Alignment: &Alignment{
					Seq:         make([]dna.Base, 0, 150),
					Path:        make([]uint32, 0, 10),
					queryStart:  0,
					targetStart: 0,
					targetEnd:   0,
					queryEnd:    0,
				},
			}
			return &pool
		},
	}
}

func NewSeedBuilder() *sync.Pool {
	return &sync.Pool{
		New: func() interface{} {
			pool := SeedMemory{
				tempSeed: Seed{},
				Hits:     make([]Seed, 0, 10000),
				Worker:   make([]Seed, 0, 10000),
			}
			return &pool
		},
	}
}

func NewMatrixPool(size int) *sync.Pool {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return &sync.Pool{
		New: func() interface{} {
			pool := Matrix{
				matrix:  m,
				trace:   trace,
				i:       0,
				j:       0,
				index:   0,
				currMax: math.MinInt64,
				route:   make([]cigar.Cigar, 0, 1),
			}
			return &pool
		},
	}
}

func NewAlignmentPool() *sync.Pool {
	return &sync.Pool{
		New: func() interface{} {
			align := Alignment{
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				queryStart:  0,
				targetStart: 0,
				targetEnd:   0,
				queryEnd:    0,
			}
			return &align
		},
	}
}

func resetDynamicScore(matrix *Matrix) {
	matrix.route = matrix.route[:0]
	matrix.currMax = math.MinInt64
}

func resetScoreKeeper(sk scoreKeeper) {
	sk.leftAlignment, sk.rightAlignment = sk.leftAlignment[:0], sk.rightAlignment[:0]
	sk.leftPath, sk.rightPath = sk.leftPath[:0], sk.rightPath[:0]
	sk.leftSeq, sk.rightSeq, sk.currSeq = sk.leftSeq[:0], sk.rightSeq[:0], sk.currSeq[:0]
}

func restartSeedHelper(helper *SeedMemory) {
	helper.currHits = helper.currHits[:0]
	helper.keyIdx, helper.keyOffset, helper.readOffset, helper.nodeOffset = 0, 0, 0, 0
	helper.nodeIdx, helper.nodePos = 0, 0
	helper.seqKey, helper.codedNodeCoord = 0, 0
	helper.leftMatches = 0
}
