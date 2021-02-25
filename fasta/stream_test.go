package fasta

import (
	"testing"
)

func TestGoReadToChan(t *testing.T) {
	for _, test := range readWriteTests {
		stream := GoReadToChan(test.filename)
		actual := make([]*Fasta, 0, len(test.data))
		for val := range stream {
			ToAppend := val
			actual = append(actual, &ToAppend)
		}
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}
