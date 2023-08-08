package obo

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ToDotTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
}{
	{InFile: "testdata/test.obo",
		OutFile:      "testdata/out.dot",
		ExpectedFile: "testdata/expected.dot"},
}

func TestToDot(t *testing.T) {
	var records map[string]*Obo
	var err error
	for _, v := range ToDotTests {
		records, _ = Read(v.InFile, true)
		ToDot(v.OutFile, records)
		if !fileio.AreEqualIgnoreOrder(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: ToDot output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var SubTreeReportTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
}{
	{InFile: "testdata/test.obo",
		OutFile:      "testdata/out.report.txt",
		ExpectedFile: "testdata/expected.report.txt"},
}

func TestSubTreeReport(t *testing.T) {
	var records map[string]*Obo
	var err error
	for _, v := range SubTreeReportTests {
		records, _ = Read(v.InFile, true)
		NumberOfDescendents(records)
		SubTreeReport(v.OutFile, records)
		if !fileio.AreEqualIgnoreOrder(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output was not as expected in SubtreeReport.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var SubTreeToDotTests = []struct {
	InFile       string
	NodeId       string
	OutFile      string
	ExpectedFile string
}{
	//{InFile: "testdata/go.obo",
	{InFile: "testdata/test.obo",
		NodeId: "GO:0000030",
		//NodeId: "GO:0042110",
		OutFile: "testdata/out.mannosyltransferaseActivity.dot",
		//OutFile:      "testdata/out.tCellActivation.dot",
		ExpectedFile: "testdata/expected.mannosyltransferaseActivity.dot",
	},
}

func TestSubtreeToDot(t *testing.T) {
	var err error
	var records map[string]*Obo
	for _, v := range SubTreeToDotTests {
		records, _ = Read(v.InFile, true)
		SubtreeToDot(v.OutFile, v.NodeId, records)
		if !fileio.AreEqualIgnoreOrder(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output not as expected in TestSubtreeToDot.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var FindTreeRootsTests = []struct {
	InFile        string
	ExpectedRoots map[string]bool
}{
	{InFile: "testdata/microTest.obo",
		ExpectedRoots: map[string]bool{"GO:1": true}},
}

func TestFindTreeRoots(t *testing.T) {
	var foundInMap bool
	var roots []*Obo
	var records map[string]*Obo
	for _, v := range FindTreeRootsTests {
		records, _ = Read(v.InFile, true)
		roots = findTreeRoots(records)
		if len(roots) != len(v.ExpectedRoots) {
			t.Errorf("Error: number of roots found: %v, did not match expected: %v.\n", len(roots), len(v.ExpectedRoots))
		}
		for i := range roots {
			if _, foundInMap = v.ExpectedRoots[roots[i].Id]; !foundInMap {
				t.Errorf("Error: root with ID: \"%v\" and Name:\"%v\" is not an expected root.\n", roots[i].Id, roots[i].Name)
			}
		}
	}
}
