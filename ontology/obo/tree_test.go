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
	var records []Obo
	var err error
	for _, v := range ToDotTests {
		records, _ = Read(v.InFile)
		ToDot(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
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
	var records []Obo
	var err error
	var termMap map[string]*Obo
	for _, v := range SubTreeReportTests {
		records, _ = Read(v.InFile)
		termMap = BuildTree(records)
		NumberOfDescendents(termMap)
		SubTreeReport(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
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
	var termMap map[string]*Obo
	var records []Obo
	for _, v := range SubTreeToDotTests {
		records, _ = Read(v.InFile)
		termMap = BuildTree(records)
		SubtreeToDot(v.OutFile, v.NodeId, termMap)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output not as expected in TestSubtreeToDot.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
