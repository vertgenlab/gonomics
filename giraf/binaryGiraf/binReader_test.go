package binaryGiraf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/giraf"
	"io"
	"reflect"
	"testing"
)

var nodeSeq []dna.Base = dna.StringToBases("ATGCGATGCG" +
	"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG" +
	"ATGCGATGCGATGCGATGCGATGCGATGCGATGCG") // 100mer of "ATGCG" repeats

func MakeTestGraph() *genomeGraph.GenomeGraph {
	var n0, n1, n2 genomeGraph.Node
	var e1, e2 genomeGraph.Edge

	n1 = genomeGraph.Node{
		Id:   1,
		Seq:  nodeSeq,
		Next: []genomeGraph.Edge{e1},
	}

	n2 = genomeGraph.Node{
		Id:   2,
		Seq:  nodeSeq,
		Prev: []genomeGraph.Edge{e2},
	}

	e1 = genomeGraph.Edge{
		Dest: &n2,
	}

	e2 = genomeGraph.Edge{
		Dest: &n1,
	}

	return &genomeGraph.GenomeGraph{
		Nodes: []genomeGraph.Node{n0, n1, n2},
	}
}

func TestRead(t *testing.T) {
	// Initialize infile
	reader := NewBinReader("testdata/test.giraf.fe")
	var err error

	// Initialize outfile
	outfile := fileio.EasyCreate("testdata/readtest.giraf")

	// Read info until EOF
	var curr giraf.Giraf
	graph := MakeTestGraph()
	for curr, err = ReadGiraf(reader, graph); err != io.EOF; curr, err = ReadGiraf(reader, graph) {
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

	err = outfile.Close()
	if err != nil {
		t.Error(err)
	}
	fileio.EasyRemove("testdata/readtest.giraf")

}

func TestReadAndWrite(t *testing.T) {
	correct := giraf.Read("testdata/test.giraf")
	CompressGiraf("testdata/test.giraf", "testdata/test.giraf.fe")

	// Initialize infile
	reader := NewBinReader("testdata/test.giraf.fe")
	var err error

	// Read info until EOF
	var curr giraf.Giraf
	graph := MakeTestGraph()
	for curr, err = ReadGiraf(reader, graph); err != io.EOF; curr, err = ReadGiraf(reader, graph) {
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
