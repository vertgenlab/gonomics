package genomeGraph

import (
	"sync"
)

type Matrix struct {
	matrix [][]int64
	trace  [][]byte
	index  int
}

func NewMatrixPool(size int) *sync.Pool {
	m := make([][]int64, size)
	t := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		t[idx] = make([]byte, size)
	}

	return &sync.Pool{
		New: func() interface{} {
			pool := Matrix{
				matrix: m,
				trace:  t,
			}
			return &pool
		},
	}
}

// Function to reset the Matrix
func (m *Matrix) Reset(rows, cols int) {
	// Check and reallocate matrix if necessary
	if cap(m.matrix) < rows {
		m.matrix = make([][]int64, rows)
		m.trace = make([][]byte, rows)
	}
	// Reset the length, but keep capacity
	m.matrix, m.trace = m.matrix[:rows], m.trace[:rows]

	for m.index = 0; m.index < rows; m.index++ {
		if cap(m.matrix[m.index]) < cols {
			m.matrix[m.index] = make([]int64, cols)
			m.trace[m.index] = make([]byte, cols)
		}
		m.matrix[m.index] = m.matrix[m.index][:cols]
		m.trace[m.index] = m.trace[m.index][:cols]

		// m.matrix[m.index][0] = 0
		// m.trace[m.index][0] = cigar.Insertion
	}
	// for m.index = 0; m.index < cols; m.index++ {
	// 	m.matrix[0][m.index] = 0
	// 	m.trace[0][m.index] = cigar.Insertion
	// }
}
