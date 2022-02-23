package main

import (
	"bytes"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"net/http"
	"os"
	"strings"
	"testing"
)

func TestVcfWebAnnotate(t *testing.T) {
	if !pingServer() {
		return
	}

	testfile := "tmp_testfile.vcf"
	vcfs, header := vcf.GoReadToChan("testdata/short.vcf")
	testout := fileio.EasyCreate(testfile)
	vcfWebAnnotate(vcfs, header, testout, 200, 2, false)
	err := testout.Close()
	exception.PanicOnErr(err)

	actualRecords, actualHeader := vcf.Read(testfile)
	expectedRecords, expectedHeader := vcf.Read(testfile)

	if len(actualRecords) != len(expectedRecords) {
		t.Error("number of records changed")
	}

	if strings.Join(actualHeader.Text, "") != strings.Join(expectedHeader.Text, "") {
		t.Error("problem writing header")
	}

	// note: I did not want to test for equality as I figure the values may change slightly as more
	// data is gathered by gnomAD and 1kgp. Therefore I am just testing to make sure annotated fields
	// exist in the output records.
	for i := range actualRecords {
		if strings.Contains(expectedRecords[i].Info, "MaxPopAF") && !strings.Contains(actualRecords[i].Info, "MaxPopAF") {
			t.Error("problem with MaxPopAF annotation")
		}

		if strings.Contains(expectedRecords[i].Info, "Consequence") && !strings.Contains(actualRecords[i].Info, "Consequence") {
			t.Error("problem with Consequence annotation")
		}

		if strings.Contains(expectedRecords[i].Info, "Gene") && !strings.Contains(actualRecords[i].Info, "Gene") {
			t.Error("problem with Gene annotation")
		}

		if strings.Contains(expectedRecords[i].Info, "Transcript") && !strings.Contains(actualRecords[i].Info, "Transcript") {
			t.Error("problem with Transcript annotation")
		}

		if strings.Contains(expectedRecords[i].Info, "ProteinEffect") && !strings.Contains(actualRecords[i].Info, "ProteinEffect") {
			t.Error("problem with ProteinEffect annotation")
		}
	}

	if !t.Failed() {
		os.Remove(testfile)
	}
}

func TestResume(t *testing.T) {
	if !pingServer() {
		return
	}

	out := fileio.EasyCreate("testdata/actual.vcf")
	vcfs, header := vcf.GoReadToChan("testdata/short.vcf")
	partial, _ := vcf.GoReadToChan("testdata/short_ann_partial.vcf")
	burnRecords(vcfs, partial)

	vcfWebAnnotate(vcfs, header, out, 1000, 2, true)
	err := out.Close()
	exception.PanicOnErr(err)

	actual, _ := vcf.Read("testdata/actual.vcf")
	expected, _ := vcf.Read("testdata/short_ann.vcf")

	if actual[len(actual)-1].Pos != expected[len(expected)-1].Pos {
		t.Error("problem with resume")
	}

	if !t.Failed() {
		os.Remove("testdata/actual.vcf")
	}
}

func pingServer() bool {
	baseUrl := "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/annotation?assembly=grch38"
	pingQuery := new(bytes.Buffer)
	pingQuery.WriteString("19:45411941:T:C")
	for i := 0; i < 200; i++ {
		pingQuery.WriteString(",19:45411941:T:C")
	}
	response, err := http.Post(baseUrl, "text/plain", pingQuery)

	if err != nil || response.StatusCode != 200 { // server not available, pass with warning
		log.Println("Could not test vcfWebAnnotate as server is unavailable")
		return false
	}
	return true
}
