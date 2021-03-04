package fasta

import (
	"testing"
)

func TestGoReadToChan(t *testing.T) {
	for _, test := range readWriteTests {
		stream := GoReadToChan(test.filename)
		actual := make([]Fasta, 0, len(test.data))
		for val := range stream {
			ToAppend := val
			actual = append(actual, ToAppend)
		}
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func sendValToChan() chan Fasta {
	answer := make(chan Fasta)

	go func(chan Fasta) {
		for i := 0; i < numRecords; i++ {
			answer <- Fasta{Name: "abcdefg"}
			//answer <- Fasta{Seq: make([]dna.Base, fastaSeqSize)}
			//answer <- Fasta{Name: "abcdefg", Seq: make([]dna.Base, fastaSeqSize)}
		}
		close(answer)
	}(answer)

	return answer
}

func sendPtrToChan() chan *Fasta {
	answer := make(chan *Fasta)

	go func(chan *Fasta) {
		for i := 0; i < numRecords; i++ {
			answer <- &Fasta{Name: "abcdefg"}
			//answer <- &Fasta{Seq: make([]dna.Base, fastaSeqSize)}
			//answer <- &Fasta{Name: "abcdefg", Seq: make([]dna.Base, fastaSeqSize)}
		}
		close(answer)
	}(answer)

	return answer
}

var numRecords int = 10000
var fastaSeqSize int = 100
var fa Fasta
var faPtr *Fasta

func BenchmarkChanCopyVal(b *testing.B) {
	for i := 0; i < b.N; i++ {
		stream := sendValToChan()
		for val := range stream {
			fa = val
		}
	}
}

func BenchmarkChanDerefPtr(b *testing.B) {
	for i := 0; i < b.N; i++ {
		stream := sendPtrToChan()
		for val := range stream {
			faPtr = val
		}
	}
}

// testing ptrs in various sized values
var kb1 [1000]byte
var kb1Ptr *[1000]byte
var kb100 [100000]byte
var kb100Ptr *[100000]byte
var mb10 [10000000]byte
var mb10Ptr *[10000000]byte

func make1kb() [1000]byte {
	var ans [1000]byte
	return ans
}

func make1kbPtr() *[1000]byte {
	var ans [1000]byte
	return &ans
}

func make100kb() [100000]byte {
	var ans [100000]byte
	return ans
}

func make100kbPtr() *[100000]byte {
	var ans [100000]byte
	return &ans
}

func make10mb() [10000000]byte {
	var ans [10000000]byte
	return ans
}

func make10mbPtr() *[10000000]byte {
	var ans [10000000]byte
	return &ans
}

func BenchmarkCopy1kb(b *testing.B) {
	for i := 0; i < b.N; i++ {
		kb1 = make1kb()
	}
}

func BenchmarkCopy1kbPtr(b *testing.B) {
	for i := 0; i < b.N; i++ {
		kb1Ptr = make1kbPtr()
	}
}

func BenchmarkCopy100kb(b *testing.B) {
	for i := 0; i < b.N; i++ {
		kb100 = make100kb()
	}
}

func BenchmarkCopy100kbPtr(b *testing.B) {
	for i := 0; i < b.N; i++ {
		kb100Ptr = make100kbPtr()
	}
}

func BenchmarkCopy10mb(b *testing.B) {
	for i := 0; i < b.N; i++ {
		mb10 = make10mb()
	}
}

func BenchmarkCopy10mbPtr(b *testing.B) {
	for i := 0; i < b.N; i++ {
		mb10Ptr = make10mbPtr()
	}
}
