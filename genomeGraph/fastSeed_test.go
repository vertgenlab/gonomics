package genomeGraph

import (
	"testing"
)

func TestCompareSeedDev(t *testing.T) {
	// Test case 1: a is equal to b
	a := &SeedDev{TargetId: 1, TargetStart: 10, Length: 20}
	b := &SeedDev{TargetId: 1, TargetStart: 10, Length: 20}
	result := CompareSeedDev(a, b)
	if result != 0 {
		t.Errorf("Error: Expected 0, got %d", result)
	}

	// Test case 2: a is less than b
	a = &SeedDev{TargetId: 1, TargetStart: 10, Length: 20}
	b = &SeedDev{TargetId: 2, TargetStart: 10, Length: 20}
	result = CompareSeedDev(a, b)
	if result != -1 {
		t.Errorf("Error: Expected -1, got %d", result)
	}

	// Test case 3: a is greater than b
	a = &SeedDev{TargetId: 2, TargetStart: 10, Length: 20}
	b = &SeedDev{TargetId: 1, TargetStart: 10, Length: 20}
	result = CompareSeedDev(a, b)
	if result != 1 {
		t.Errorf("Expected 1, got %d", result)
	}

	// Test case 4: a is equal to b in some fields but not all
	a = &SeedDev{TargetId: 1, TargetStart: 10, Length: 20}
	b = &SeedDev{TargetId: 1, TargetStart: 10, Length: 30}
	result = CompareSeedDev(a, b)
	if result != -1 {
		t.Errorf("Error: Expected -1, got %d", result)
	}

	// Test case 5: a and b have different values in all fields
	a = &SeedDev{TargetId: 1, TargetStart: 10, Length: 20}
	b = &SeedDev{TargetId: 2, TargetStart: 20, Length: 30}
	result = CompareSeedDev(a, b)
	if result != -1 {
		t.Errorf("Error: Expected -1, got %d", result)
	}
}
