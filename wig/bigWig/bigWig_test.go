package bigWig

import (
	"testing"
)

var readHeaderTests = []struct {
	InFile               string
	ExpectedBbiHeader    BbiHeader
	ExpectedZoomHeaders  []ZoomHeader
	ExpectedTotalSummary TotalSummaryBlock
}{
	{InFile: "testdata/test.bw",
		ExpectedBbiHeader: BbiHeader{
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
		},
		ExpectedZoomHeaders: []ZoomHeader{
			{ReductionLevel: 66, Reserved: 0, DataOffset: 6457, IndexOffset: 6492},
			{ReductionLevel: 264, Reserved: 0, DataOffset: 12696, IndexOffset: 12731},
		},
		ExpectedTotalSummary: TotalSummaryBlock{
			BasesCovered: 15,
			MinVal:       6,
			MaxVal:       47,
			SumData:      208,
			SumSquares:   4144,
		},
	},
}

func TestRead(t *testing.T) {
	var answer BigWig
	var currZoomHeader int
	for _, v := range readHeaderTests {
		answer = Read(v.InFile)
		testBbiHeaders(answer.BbiHeader, v.ExpectedBbiHeader, t)
		for currZoomHeader = 0; currZoomHeader < len(v.ExpectedZoomHeaders); currZoomHeader++ {
			testZoomHeaders(answer.ZoomHeaders[currZoomHeader], v.ExpectedZoomHeaders[currZoomHeader], t)
		}
		testTotalSummary(answer.TotalSummaryBlock, v.ExpectedTotalSummary, t)
	}
}

func testBbiHeaders(a BbiHeader, b BbiHeader, t *testing.T) {
	if a.Magic != b.Magic {
		t.Errorf("Error: header magic not as expected.\n")
	}
	if a.Version != b.Version {
		t.Errorf("Error: header version not as expected.\n")
	}
	if a.ZoomLevels != b.ZoomLevels {
		t.Errorf("Error: header zoom levels not as expected.\n")
	}
	if a.ChromosomeTreeOffset != b.ChromosomeTreeOffset {
		t.Errorf("Error: header chromosome tree offset not  as expected.\n")
	}
	if a.FullDataOffset != b.FullDataOffset {
		t.Errorf("Error: header full data offset not as expected.\n")
	}
	if a.FieldCount != b.FieldCount {
		t.Errorf("Error: header field count was not as expected.\n")
	}
	if a.DefinedFieldCount != b.DefinedFieldCount {
		t.Errorf("Error: header defined field count was not as expected.\n")
	}
	if a.AutoSqlOffset != b.AutoSqlOffset {
		t.Errorf("Error: header autoSqlOffset field was not as expected.\n")
	}
	if a.TotalSummaryOffset != b.TotalSummaryOffset {
		t.Errorf("Error: header TotalSummaryOffset field was not as expected.\n")
	}
	if a.UncompressBufferSize != b.UncompressBufferSize {
		t.Errorf("Error: header UncompressBufferSize field was not as expected.\n")
	}
	if a.ExtensionOffset != b.ExtensionOffset {
		t.Errorf("Error: heaader ExtensionOffset field was not as expected.\n")
	}
}

func testZoomHeaders(a ZoomHeader, b ZoomHeader, t *testing.T) {
	if a.ReductionLevel != b.ReductionLevel {
		t.Errorf("Error: ZoomHeader ReductionLevel field was not as expected.\n")
	}
	if a.Reserved != b.Reserved {
		t.Errorf("Error: ZoomHeader Reserved field was not as expected.\n")
	}
	if a.DataOffset != b.DataOffset {
		t.Errorf("Error: ZoomHeader DataOffset field was not as expected.\n")
	}
	if a.IndexOffset != b.IndexOffset {
		t.Errorf("Error: ZoomHeader IndexOffset field was not as expected.\n")
	}
}

func testTotalSummary(a TotalSummaryBlock, b TotalSummaryBlock, t *testing.T) {
	if a.BasesCovered != b.BasesCovered {
		t.Errorf("Error: TotalSummaryBlock BasesCovered field was not as expected.\n")
	}
	if a.MinVal != b.MinVal {
		t.Errorf("Error: TotalSummaryBlock MinVal field was not as expected.\n")
	}
	if a.MaxVal != b.MaxVal {
		t.Errorf("Error: TotalSummaryBlock MaxVal field was not as expected.\n")
	}
	if a.SumData != b.SumData {
		t.Errorf("Error: TotalSummaryBlock SumData field was not as expected.\n")
	}
	if a.SumSquares != b.SumSquares {
		t.Errorf("Error: TotalSummaryBlock SumSquares field was not as expected.\n")
	}
}
