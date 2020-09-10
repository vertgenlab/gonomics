package giraf

import (
	"bytes"
	"github.com/biogo/hts/bgzf"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"io"
	"strings"
)

type BinReader struct {
	bg  *bgzf.Reader
	buf bytes.Buffer
}

func NewBinReader(file io.Reader) *BinReader {
	reader, err := bgzf.NewReader(file, 1) //TODO: Play with different levels of concurrency
	common.ExitIfError(err)
	return &BinReader{
		bg: reader,
	}
}

func DecompressGiraf(filename string, graph *simpleGraph.SimpleGraph) {
	// Initialize infile
	infile := fileio.EasyOpen(filename)
	defer infile.Close()
	reader := NewBinReader(infile.BuffReader)
	var err error

	// Initialize oufile
	outfile := fileio.EasyCreate(strings.TrimSuffix(filename, ".fe"))
	defer outfile.Close()

	// Read info until EOF
	var curr Giraf
	for curr, err = reader.Read(graph); err != io.EOF; curr, err = reader.Read(graph) {
		common.ExitIfError(err)
		WriteGirafToFileHandle(outfile, &curr)
	}

	// Close reader
	err = reader.bg.Close()
	common.ExitIfError(err)
}

func (br *BinReader) Read(g *simpleGraph.SimpleGraph) (Giraf, error) {
	var answer Giraf
	var err error

	return answer, err
}
