package binaryGiraf

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"io"
	"testing"
)

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
	for curr, err = reader.Read(&simpleGraph.SimpleGraph{}); err != io.EOF; curr, err = reader.Read(&simpleGraph.SimpleGraph{}) {
		common.ExitIfError(err)
		giraf.WriteGirafToFileHandle(outfile, &curr)
	}

	// Close reader
	err = reader.bg.Close()
	common.ExitIfError(err)
}
