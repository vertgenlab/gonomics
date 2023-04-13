package main

import (
	"encoding/csv"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var VcfFormatTests = []struct {
	InFile        string
	OutFile       string
	ExpectedFile  string
	EnsemblToUCSC bool
	UCSCToEnsembl bool
	FixVcfRecords bool
	Ref           string
	ClearInfo     bool
}{
	{"testdata/test.UCSC.vcf", "testdata/tmp.UCSCtoEnsembl.vcf", "testdata/test.Ensembl.vcf", false, true, false, "", false},
	{"testdata/test.Ensembl.vcf", "testdata/tmp.EnsemblToUCSC.vcf", "testdata/test.UCSC.vcf", true, false, false, "", false},
	{"testdata/test.UCSC.vcf", "testdata/tmp.UCSCnoInfo.vcf", "testdata/expected.noInfo.vcf", false, false, false, "", true},
	{"testdata/test.broken.vcf", "testdata/tmp.fixed.vcf", "testdata/expected.fixed.vcf", false, false, true, "testdata/test.fa", false},
}

func TestVcfFormat(t *testing.T) {
	var err error
	for _, v := range VcfFormatTests {
		vcfFormat(v.InFile, v.OutFile, v.EnsemblToUCSC, v.UCSCToEnsembl, v.FixVcfRecords, v.Ref, v.ClearInfo, false)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in VcfFormat. Output does not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

// TODO: better table tests.
func TestVcfTable(t *testing.T) {
	infile := "testdata/test_table.vcf"
	outfile := "testdata/actual_table.csv"
	expected := "testdata/table_expected.csv"
	vcfFormat(infile, outfile, false, false, false, "", false, true)

	actualFile := fileio.EasyOpen(outfile)
	expectedFile := fileio.EasyOpen(expected)
	actualReader := csv.NewReader(actualFile)
	expectedReader := csv.NewReader(expectedFile)

	var actualData, expectedData [][]string
	var err error
	expectedData, err = expectedReader.ReadAll()
	exception.PanicOnErr(err)
	actualData, err = actualReader.ReadAll()
	if err != nil {
		t.Error("problem writing csv file")
	}

	var actualOrdering []int = make([]int, len(expectedData[0]))

	// get order to deal with map randomness
	for i := range expectedData[0] {
		for j := range actualData[0] {
			if expectedData[0][i] == actualData[0][j] {
				actualOrdering[j] = i
			}
		}
	}

	// reorder actual
	for i := range actualData {
		newLine := make([]string, len(actualData[i]))
		for j := range actualData[i] {
			newLine[actualOrdering[j]] = actualData[i][j]
		}
		actualData[i] = newLine
	}

	if len(expectedData) != len(actualData) {
		t.Error("problem writing csv file")
	}

	for i := range expectedData {
		if strings.Join(expectedData[i], "") != strings.Join(actualData[i], "") {
			t.Error("problem with vcf formatting to csv")
		}
	}

	err = actualFile.Close()
	exception.PanicOnErr(err)
	err = expectedFile.Close()
	exception.PanicOnErr(err)

	if !t.Failed() {
		err = os.Remove(outfile)
		exception.PanicOnErr(err)
	}
}
