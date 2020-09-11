package binaryGiraf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"io"
	"reflect"
	"testing"
)

var nodeSeq []dna.Base = dna.StringToBases("ATGCGATGCG" +
	"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG" +
	"ATGCGATGCGATGCGATGCGATGCGATGCGATGCG") // 100mer of "ATGCG" repeats

func MakeTestGraph() *simpleGraph.SimpleGraph {
	var n0, n1, n2 simpleGraph.Node
	var e1, e2 simpleGraph.Edge

	n1 = simpleGraph.Node{
		Id:   1,
		Seq:  nodeSeq,
		Next: []*simpleGraph.Edge{&e1},
	}

	n2 = simpleGraph.Node{
		Id:   2,
		Seq:  nodeSeq,
		Prev: []*simpleGraph.Edge{&e2},
	}

	e1 = simpleGraph.Edge{
		Dest: &n2,
	}

	e2 = simpleGraph.Edge{
		Dest: &n1,
	}

	return &simpleGraph.SimpleGraph{
		Nodes: []*simpleGraph.Node{&n0, &n1, &n2},
	}
}

func TestRead(t *testing.T) {
	// Initialize infile
	infile := fileio.EasyOpen("testdata/test.giraf.fe")
	defer infile.Close()
	reader := NewBinReader(infile.BuffReader)
	var err error

	// Initialize outfile
	outfile := fileio.EasyCreate("testdata/readtest.giraf")
	defer outfile.Close()

	// Read info until EOF
	var curr giraf.Giraf
	graph := MakeTestGraph()
	for curr, err = reader.Read(graph); err != io.EOF; curr, err = reader.Read(graph) {
		if err != nil {
			t.Error(err)
		}
		giraf.WriteGirafToFileHandle(outfile, &curr)
	}

	// Close reader
	err = reader.bg.Close()
	if err != nil {
		t.Error(err)
	}
}

func TestReadAndWrite(t *testing.T) {
	correct := giraf.Read("testdata/test.giraf")
	CompressGiraf("testdata/test.giraf")

	// Initialize infile
	infile := fileio.EasyOpen("testdata/test.giraf.fe")
	defer infile.Close()
	reader := NewBinReader(infile.BuffReader)
	var err error

	// Read info until EOF
	var curr giraf.Giraf
	graph := MakeTestGraph()
	for curr, err = reader.Read(graph); err != io.EOF; curr, err = reader.Read(graph) {
		if !reflect.DeepEqual(*(correct[0]), curr) {
			t.Error("Problem with writing and reading binary giraf files")
		}
		if err != nil {
			t.Error(err)
		}
	}

	// Close reader
	err = reader.bg.Close()
	if err != nil {
		t.Error(err)
	}
}
