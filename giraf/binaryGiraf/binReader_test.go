package binaryGiraf

import (
	"github.com/vertgenlab/gonomics/simpleGraph"
	"testing"
)

func TestRead(t *testing.T) {
	DecompressGiraf("testdata/test.giraf.fe", &simpleGraph.SimpleGraph{})
}