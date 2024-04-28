package bigWig

import (
	"testing"
)

var readHeaderTests = []struct {
	InFile                  string
	ExpectedBbiHeader       BbiHeader
	ExpectedZoomHeaders     []ZoomHeader
	ExpectedTotalSummary    TotalSummaryBlock
	ExpectedChromTreeHeader ChromTreeHeader
	ExpectedChromTreeNodes  []ChromTreeNode
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
		ExpectedChromTreeHeader: ChromTreeHeader{
			Magic:     2026540177,
			BlockSize: 1,
			KeySize:   4,
			ValSize:   8,
			ItemCount: 1,
			Reserved:  0,
		},
		ExpectedChromTreeNodes: []ChromTreeNode{
			{IsLeaf: true,
				Reserved: 0,
				Count:    1,
				Items: []ChromTreeItem{
					{Key: []byte{99, 104, 114, 49}, // this spells "chr1" in ASCII.
						ChromId:     0,
						ChromSize:   20,
						ChildOffset: 0,
					},
				},
			},
		},
	},
	{InFile: "testdata/wholeGenome.bw",
		ExpectedBbiHeader: BbiHeader{
			Magic:                bigWigMagic,
			Version:              4,
			ZoomLevels:           2,
			ChromosomeTreeOffset: 152,
			FullDataOffset:       212,
			FullIndexOffset:      322,
			FieldCount:           0, // this value must always be zero for valid bigWigs.
			DefinedFieldCount:    0, // this value must be zero for valid bigWigs.
			AutoSqlOffset:        0, //this value is unused in bigWigs as well.
			TotalSummaryOffset:   112,
			UncompressBufferSize: 32768,
			ExtensionOffset:      0,
		},
		ExpectedZoomHeaders: []ZoomHeader{
			{ReductionLevel: 30, Reserved: 0, DataOffset: 6542, IndexOffset: 6613},
			{ReductionLevel: 120, Reserved: 0, DataOffset: 12817, IndexOffset: 12870},
		},
		ExpectedTotalSummary: TotalSummaryBlock{
			BasesCovered: 17,
			MinVal:       2,
			MaxVal:       12,
			SumData:      104,
			SumSquares:   830,
		},
		ExpectedChromTreeHeader: ChromTreeHeader{
			Magic:     2026540177,
			BlockSize: 2,
			KeySize:   4,
			ValSize:   8,
			ItemCount: 2,
			Reserved:  0,
		},
		ExpectedChromTreeNodes: []ChromTreeNode{
			{IsLeaf: true,
				Reserved: 0,
				Count:    2,
				Items: []ChromTreeItem{
					{Key: []byte{99, 104, 114, 65}, // this spells "chrA" in ASCII.
						ChromId:     0,
						ChromSize:   50,
						ChildOffset: 0,
					},
					{Key: []byte{99, 104, 114, 66}, // this spells "chrB" in ASCII.
						ChromId:     1,
						ChromSize:   20,
						ChildOffset: 0,
					},
				},
			},
			{IsLeaf: true,
				Reserved: 0,
				Count:    0,
			},
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
		testChromTreeHeader(answer.ChromTreeHeader, v.ExpectedChromTreeHeader, t)
		for currNodeIdx := range answer.ChromTreeNodes {
			testChromTreeNode(answer.ChromTreeNodes[currNodeIdx], v.ExpectedChromTreeNodes[currNodeIdx], t)
		}
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
	if a.FullIndexOffset != b.FullIndexOffset {
		t.Errorf("Error: header full index offset not as expected.\n")
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

func testChromTreeHeader(a ChromTreeHeader, b ChromTreeHeader, t *testing.T) {
	if a.Magic != b.Magic {
		t.Errorf("Error: ChromTreeHeader Magic field was not as expected.\n")
	}
	if a.BlockSize != b.BlockSize {
		t.Errorf("Error: ChromTreeHeader BlockSize field was not as expected.\n")
	}
	if a.KeySize != b.KeySize {
		t.Errorf("Error: ChromTreeHeader KeySize field was not as expected.\n")
	}
	if a.ValSize != b.ValSize {
		t.Errorf("Error: ChromTreeHeader ValSize field was not as expected.\n")
	}
	if a.ItemCount != b.ItemCount {
		t.Errorf("Error: ChromTreeHeader ItemCount field was not as expected.\n")
	}
	if a.Reserved != b.Reserved {
		t.Errorf("Error: ChromTreeHeader Reserved field was not as expected.\n")
	}
}

func testChromTreeNode(a ChromTreeNode, b ChromTreeNode, t *testing.T) {
	if a.IsLeaf != b.IsLeaf {
		t.Errorf("Error: ChromTreeNode IsLeaf field is not as expected.\n")
	}
	if a.Reserved != b.Reserved {
		t.Errorf("Error: ChromTreeNode Reserved field is not as expected.\n")
	}
	if a.Count != b.Count {
		t.Errorf("Error: ChromTreeNode Count field is not as expected.\n")
	}
	if len(a.Items) != len(b.Items) {
		t.Errorf("Error: ChromTreeNode does not have the expected number of items.\n")
	}
	for currItemIdx := range a.Items {
		testChromTreeItem(a.Items[currItemIdx], b.Items[currItemIdx], t)
	}
}

func testChromTreeItem(a ChromTreeItem, b ChromTreeItem, t *testing.T) {
	if len(a.Key) != len(b.Key) {
		t.Errorf("Error: ChromTreeItem Key length is not as expected.\n")
	}
	for currKeyPos := range a.Key {
		if a.Key[currKeyPos] != b.Key[currKeyPos] {
			t.Errorf("Error: CurrTreeItem Key entry is not as expected.\n")
		}
		if a.ChromId != b.ChromId {
			t.Errorf("Error: CurrTreeItem ChromId entry is not as expected.\n")
		}
		if a.ChromSize != b.ChromSize {
			t.Errorf("Error: CurrTreeItem ChromSize entry is not as expected.\n")
		}
		if a.ChildOffset != b.ChildOffset {
			t.Errorf("Error: CurrTreeItem ChildOffset entry is not as expected.\n")
		}
	}
}
