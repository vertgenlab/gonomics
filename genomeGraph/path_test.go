package genomeGraph

import (
	"testing"
)

func TestCatPaths(t *testing.T) {
	cases := []struct {
		name      string
		currPaths []uint32
		newPaths  []uint32
		expected  []uint32
	}{
		{
			name:      "Empty Current Paths and New Paths",
			currPaths: []uint32{},
			newPaths:  []uint32{},
			expected:  []uint32{},
		},
		{
			name:      "Empty Current Paths",
			currPaths: []uint32{},
			newPaths:  []uint32{1, 2, 3},
			expected:  []uint32{1, 2, 3},
		},
		{
			name:      "Empty New Paths",
			currPaths: []uint32{1, 2, 3},
			newPaths:  []uint32{},
			expected:  []uint32{1, 2, 3},
		},
		{
			name:      "Non-empty Current Paths and New Paths",
			currPaths: []uint32{1, 2, 3},
			newPaths:  []uint32{4, 5, 6},
			expected:  []uint32{1, 2, 3, 4, 5, 6},
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			result := CatPaths(c.currPaths, c.newPaths)
			if len(result) != len(c.expected) {
				t.Errorf("Expected %v, got %v", c.expected, result)
			} else {
				for i := 0; i < len(c.expected); i++ {
					if result[i] != c.expected[i] {
						t.Errorf("Paths are not equal expected %v, got %v", c.expected, result)
					}
				}
			}
		})
	}
}
func TestAddPath(t *testing.T) {
	cases := []struct {
		name     string
		allPaths []uint32
		newPath  uint32
		expected []uint32
	}{
		{
			name:     "Empty All Paths",
			allPaths: []uint32{},
			newPath:  1,
			expected: []uint32{1},
		},
		{
			name:     "All Paths Contains New Path",
			allPaths: []uint32{1, 2, 3},
			newPath:  3,
			expected: []uint32{1, 2, 3},
		},
		{
			name:     "All Paths Does Not Contain New Path",
			allPaths: []uint32{1, 2, 3},
			newPath:  4,
			expected: []uint32{1, 2, 3, 4},
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			result := AddPath(c.allPaths, c.newPath)
			if len(result) != len(c.expected) {
				t.Errorf("Expected %v, got %v", c.expected, result)
			} else {
				for i := 0; i < len(c.expected); i++ {
					if result[i] != c.expected[i] {
						t.Errorf("Paths are not equal expected %v, got %v", c.expected, result)
					}
				}
			}
		})
	}
}
