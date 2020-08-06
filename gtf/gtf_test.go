package gtf

import (
	"os"
	"testing"
)

var c1 CDS = CDS{Start: 173477398, End: 173477492, Score: -1, Frame: 0}
var f FiveUTR = FiveUTR{Start: 173477335, End: 173477397, Score: -1}

//var e1 Exon = Exon{Start: 173477335, End: 173477492, Score: -1, ExonNumber: "\"1\"", ExonID: "\"NM_004905.1\"", Cds: &c1, FiveUtr: &f}

var c2 CDS = CDS{Start: 173487735, End: 173487860, Score: -1, Frame: 0}
var t ThreeUTR = ThreeUTR{Start: 173487864, End: 173488815, Score: -1}

//var e2 Exon = Exon{Start: 173487735, End: 173488815, Score: -1, ExonNumber: "\"5\"", ExonID: "\"NM_004905.5\"", Cds: &c2, ThreeUtr: &t}
//var t1 = Transcript{Chr: "chr1", Source: "refGene", Start: 173477335, End: 173488815, Score: -1, Strand: true, TranscriptID: "\"NM_004905\"", Exons: []*Exon{&e1, &e2}}
//var g1 Gene = Gene{GeneID: "\"PRDX6\"", GeneName: "\"PRDX6\"", Transcripts: []*Transcript{&t1}}

var e1Internal Exon = Exon{Start: 173477335, End: 173477492, Score: -1, ExonNumber: "1", ExonID: "NM_004905.1", Cds: &c1, FiveUtr: &f}
var e2Internal Exon = Exon{Start: 173487735, End: 173488815, Score: -1, ExonNumber: "5", ExonID: "NM_004905.5", Cds: &c2, ThreeUtr: &t}
var t1Internal = Transcript{Chr: "chr1", Source: "refGene", Start: 173477335, End: 173488815, Score: -1, Strand: true, TranscriptID: "NM_004905", Exons: []*Exon{&e1Internal, &e2Internal}}
var g1Internal Gene = Gene{GeneID: "PRDX6", GeneName: "PRDX6", Transcripts: []*Transcript{&t1Internal}}

//, &g2}
/*
var readWriteTests = []struct {
	filename string
	data     map[string]*Gene
}{
	{"testdata/gtfFileTest.gtf", g},
}*/

func TestRead(t *testing.T) {
	var g map[string]*Gene = make(map[string]*Gene)
	g["PRDX6"] = &g1Internal
	actual := Read("testdata/gtfFileTest.gtf")
	if !AllAreEqual(actual, g) {
		t.Errorf("The file was not read correctly.")
	}
	/*
		for _, test := range readWriteTests {
			actual := Read(test.filename)
			if !AllAreEqual(actual, test.data) {
				t.Errorf("The %s file was not read correctly.", test.filename)
			}
		}*/
}

func TestWriteAndRead(t *testing.T) {
	var g map[string]*Gene = make(map[string]*Gene)
	g["PRDX6"] = &g1Internal
	tempFile := "testdata/gtfFileTest.gtf.tmp"
	Write(tempFile, g)
	whatIRead := Read(tempFile)
	if !AllAreEqual(g, whatIRead) {
		t.Errorf("The %s file was not written and read correctly.", "testdata/gtfFileTest.gtf")
	}
	err := os.Remove(tempFile)
	if err != nil {
		t.Errorf("Deleting temp file %s gave an error.", tempFile)
	}
}

/*
func TestWriteAndRead(t *testing.T) {
	var actual map[string]*Gene
	var g map[string]*Gene
	g["\"PRDX6\""] = &g1
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		Write(tempFile, test.data)
		actual = Read(tempFile)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not written and read correctly.", test.filename)
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}*/
