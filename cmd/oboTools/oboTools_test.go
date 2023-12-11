package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var OboToolsMappingTests = []struct {
	InFile       string
	OutFile      string
	Force        bool
	ExpectedFile string
}{
	{InFile: "../../ontology/obo/testdata/test.obo",
		OutFile:      "testdata/out.mapping.txt",
		Force:        true,
		ExpectedFile: "testdata/expected.mapping.txt",
	},
}

func TestOboToolsMapping(t *testing.T) {
	var err error
	var s MappingSettings
	for _, v := range OboToolsMappingTests {
		s = MappingSettings{
			InFile:  v.InFile,
			OutFile: v.OutFile,
			Force:   v.Force,
		}
		OboToolsMapping(s)
		if !fileio.AreEqualIgnoreOrder(v.OutFile, v.ExpectedFile) { //because this is written from a map, order is not deterministic
			t.Errorf("Error: OboToolsMapping outfile was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
