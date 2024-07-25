package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
)

func TestNewMatrixPoolAndReset(t *testing.T) {
	poolSize := 5 // Example size for the matrix pool

	// Create the matrix pool
	pool := NewMatrixPool(poolSize)

	// Get a matrix from the pool
	matrix := pool.Get().(*Matrix)

	// Check initial dimensions
	if len(matrix.matrix) != poolSize || len(matrix.trace) != poolSize {
		t.Errorf("Expected initial dimensions %d x %d, got %d x %d", poolSize, poolSize, len(matrix.matrix), len(matrix.trace))
	}

	// Check if all elements are initialized to zero
	for i := 0; i < poolSize; i++ {
		for j := 0; j < poolSize; j++ {
			if matrix.matrix[i][j] != 0 {
				t.Errorf("Matrix element [%d][%d] not initialized to 0", i, j)
			}
			if matrix.trace[i][j] != 0 {
				t.Errorf("Trace element [%d][%d] not initialized to 0", i, j)
			}
		}
	}

	// Reset the matrix to different dimensions (smaller)
	newRows, newCols := 3, 4
	matrix.Reset(newRows, newCols)

	// Check new dimensions
	if len(matrix.matrix) != newRows || len(matrix.trace) != newRows {
		t.Errorf("Expected dimensions after reset: %d x %d, got %d x %d", newRows, newCols, len(matrix.matrix), len(matrix.trace))
	}
	for _, row := range matrix.matrix {
		if len(row) != newCols {
			t.Errorf("Expected row length %d, got %d", newCols, len(row))
		}
	}
	for _, row := range matrix.trace {
		if len(row) != newCols {
			t.Errorf("Expected row length %d, got %d", newCols, len(row))
		}
	}

	// Check if the first row and column are set correctly
	for i := 0; i < newRows; i++ {
		if matrix.matrix[i][0] != 0 {
			t.Errorf("Matrix element [%d][0] not set to 0 after reset", i)
		}
		if matrix.trace[i][0] != cigar.Insertion {
			t.Errorf("Trace element [%d][0] not set to cigar.Insertion after reset", i)
		}
	}

	for j := 0; j < newCols; j++ {
		if matrix.matrix[0][j] != 0 {
			t.Errorf("Matrix element [0][%d] not set to 0 after reset", j)
		}
		if matrix.trace[0][j] != cigar.Insertion {
			t.Errorf("Trace element [0][%d] not set to cigar.Insertion after reset", j)
		}
	}
}
