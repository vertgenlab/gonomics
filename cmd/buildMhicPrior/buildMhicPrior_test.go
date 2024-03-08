package main

import "testing"

func TestBisectLeft(t *testing.T) {
	var slc []float64 = []float64{1, 3.3, 6.1, 7, 9.9, 10, 11.45, 15, 20.9}
	if bisectLeft(slc, 6.9) != 3 {
		t.Errorf("error in bisect left. Expected 3, got %d", bisectLeft(slc, 6.9))
	} else if bisectLeft(slc, 10) != 5 {
		t.Errorf("error in bisect left. Expected 5, got %d", bisectLeft(slc, 10))
	} else if bisectLeft(slc, 100) != len(slc) {
		t.Errorf("error in bisect left. Expected %d, got %d", len(slc), bisectLeft(slc, 10))
	}
}
