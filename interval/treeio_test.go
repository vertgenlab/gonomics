package interval

import (
	"bytes"
	"encoding/gob"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"testing"
)

func TestReadWriteTree(t *testing.T) {
	tree := makeTestTree()
	encodedTree := new(bytes.Buffer)
	err := WriteTree(encodedTree, tree)
	if err != nil {
		t.Errorf("problem writing tree: %s", err)
	}
	readTree, err := ReadTree(encodedTree)
	if err != nil {
		t.Errorf("problem reading tree: %s", err)
	}

	queryIntervals := generateIntervals(1000, 0, 10000, 4)
	var actualAns, testAns []Interval
	for i := range queryIntervals {
		actualAns = Query(tree, queryIntervals[i], "any")
		testAns = Query(readTree, queryIntervals[i], "any")
		if len(actualAns) != len(testAns) {
			t.Error("problem with tree encoding/decoding")
		}
		for j := range actualAns {
			if !bed.Equal(*actualAns[j].(*bed.Bed), *testAns[j].(*bed.Bed)) {
				t.Error("problem with tree encoding/decoding")
			}
		}
	}
}

func makeTestTree() Tree {
	var intervals []Interval
	intervals = append(intervals, generateIntervals(10, 0, 10000, 1)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 2)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 3)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 4)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 5)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 6)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 7)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 8)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 9)...)
	intervals = append(intervals, generateIntervals(10, 0, 10000, 10)...)
	return BuildTree(intervals)
}

func BenchmarkRegister(b *testing.B) {
	tree := makeTestTree()
	for i := 0; i < b.N; i++ {
		for _, val := range tree {
			for i := range val.Data {
				gob.Register(val.Data[i])
			}
		}
	}
}

func BenchmarkRegisterSmall(b *testing.B) {
	tree := makeTestTree()
	for i := 0; i < b.N; i++ {
		for _, val := range tree {
			for i := range val.Data {
				gob.Register(val.Data[i])
			}
			break
		}
	}
}

func BenchmarkRegisterSingle(b *testing.B) {
	tree := makeTestTree()
	for i := 0; i < b.N; i++ {
		for _, val := range tree {
			gob.Register(val.Data[0])
			break
		}
	}
}

func BenchmarkRegisterDouble(b *testing.B) {
	tree := makeTestTree()
	for i := 0; i < b.N; i++ {
		for _, val := range tree {
			gob.Register(val.Data[0])
			gob.Register(val.Data[1])
			break
		}
	}
}

func BenchmarkWriteTree(b *testing.B) {
	tree := makeTestTree()
	var err error
	for i := 0; i < b.N; i++ {
		err = WriteTree(io.Discard, tree)
		exception.PanicOnErr(err)
	}
}
