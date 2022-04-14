package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SlurmGeneralizedCheckTests = []struct {
	inputFancyFile					string
	actualOutputParseTheInput		string
	expectedOutputParseTheInput		string
} {
	{"testdata/inputFancyFile.txt", "testdata/actualOutputParseTheInput.txt", "testdata/expectedOutputParseTheInput.txt"},
}

func TestParseTheInput(t *testing.T) {
	for _,v := range SlurmGeneralizedCheckTests {
		parsed := parseTheInput(v.inputFancyFile)
		begin := parsed[0].begin
		out := parsed[0].outToCheck
		check := parsed[0].checkType
		end := parsed[0].end
		actual, err := os.Create(v.actualOutputParseTheInput)
		if err != nil {
			fmt.Println(err)
			return
		}
		fmt.Fprintf(actual,"begin: %s \n out: %s \n check: %s \n end: %s \n", begin, out, check, end)
		if !fileio.AreEqual(v.actualOutputParseTheInput, v.expectedOutputParseTheInput) {
			t.Errorf("Error in slurmGeneralizedCheck, expected %s, actual %s", v.expectedOutputParseTheInput, v.actualOutputParseTheInput)
		} else {
			err = os.Remove(v.actualOutputParseTheInput)
			exception.PanicOnErr(err)
		}

	}

}