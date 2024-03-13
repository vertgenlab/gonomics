package bigWig

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var readHeaderTests = []struct {
	InFile         string
	ExpectedHeader Header
}{
	{InFile: "testdata/test.bw",
		ExpectedHeader: Header{
			Magic:                bigWigMagic,
			Version:              4,
			ZoomLevels:           2,
			ChromosomeTreeOffset: 152,
			FullDataOffset:       200,
			FullIndexOffset:      253,
			FieldCount:           0, // this value must always be zero for valid bigWigs.
			DefinedFieldCount:    0, // this value must be zero for valid bigWigs.
			AutoSqlOffset:        0, //this value is unused in bigWigs as well.
			TotalSummaryOffset:   112,
			UncompressBufferSize: 32768,
			ExtensionOffset:      0,
		}},
}

func TestReadHeader(t *testing.T) {
	var file *fileio.EasyReader
	var header Header
	var err error
	for _, v := range readHeaderTests {
		file = fileio.EasyOpen(v.InFile)
		header = readHeader(file)
		err = file.Close()
		exception.PanicOnErr(err)
		if header.Magic != v.ExpectedHeader.Magic {
			t.Errorf("Error: header magic not as expected.\n")
		}
		if header.Version != v.ExpectedHeader.Version {
			t.Errorf("Error: header version not as expected.\n")
		}
		if header.ZoomLevels != v.ExpectedHeader.ZoomLevels {
			t.Errorf("Error: header zoom levels not as expected.\n")
		}
		if header.ChromosomeTreeOffset != v.ExpectedHeader.ChromosomeTreeOffset {
			t.Errorf("Error: header chromosome tree offset not  as expected.\n")
		}
		if header.FullDataOffset != v.ExpectedHeader.FullDataOffset {
			t.Errorf("Error: header full data offset not as expected.\n")
		}
		if header.FieldCount != v.ExpectedHeader.FieldCount {
			t.Errorf("Error: header field count was not as expected.\n")
		}
		if header.DefinedFieldCount != v.ExpectedHeader.DefinedFieldCount {
			t.Errorf("Error: header defined field count was not as expected.\n")
		}
		if header.AutoSqlOffset != v.ExpectedHeader.AutoSqlOffset {
			t.Errorf("Error: header autoSqlOffset field was not as expected.\n")
		}
		if header.TotalSummaryOffset != v.ExpectedHeader.TotalSummaryOffset {
			t.Errorf("Error: header TotalSummaryOffset field was not as expected.\n")
		}
		if header.UncompressBufferSize != v.ExpectedHeader.UncompressBufferSize {
			t.Errorf("Error: header UncompressBufferSize field was not as expected.\n")
		}
		if header.ExtensionOffset != v.ExpectedHeader.ExtensionOffset {
			t.Errorf("Error: heaader ExtensionOffset field was not as expected.\n")
		}
	}

}
