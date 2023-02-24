package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"math/rand"
	"os"
	"testing"
)

var samConsensusTests = []struct {
	inFile                     string
	refFile                    string
	outFile                    string
	multiFaDir                 string
	substitutionsOnly          bool
	insertionThreshold         float64
	tName                      string
	qName                      string
	outFile_expected           string
	multiFaOutAndExpectedFiles map[string]string //maps output file names to the expected file names
}{
	{inFile: "testdata/test.sam",
		refFile:            "testdata/test.ref.fa",
		outFile:            "testdata/tmpOutFile.fa",
		multiFaDir:         "",
		substitutionsOnly:  true,
		insertionThreshold: 0.9,
		tName:              "target",
		qName:              "query",
		outFile_expected:   "testdata/test.out.fa"},
	{inFile: "testdata/test.sam",
		refFile:                    "testdata/test.ref.fa",
		outFile:                    "testdata/tmpOutFile.indel.fa",
		multiFaDir:                 "testdata/multiFa",
		substitutionsOnly:          false,
		insertionThreshold:         0.9,
		tName:                      "target",
		qName:                      "query",
		outFile_expected:           "testdata/test.out.indel.fa",
		multiFaOutAndExpectedFiles: map[string]string{"testdata/multiFa/chr1.fa": "testdata/multiFa/expected.chr1.fa", "testdata/multiFa/chr2.fa": "testdata/multiFa/expected.chr2.fa"},
	},
}

func TestSamConsensus(t *testing.T) {
	rand.Seed(1)
	var err error
	var s Settings
	for _, v := range samConsensusTests {
		s = Settings{
			SamFileName:        v.inFile,
			RefFile:            v.refFile,
			OutFile:            v.outFile,
			MultiFaDir:         v.multiFaDir,
			SubstitutionsOnly:  v.substitutionsOnly,
			InsertionThreshold: v.insertionThreshold,
			tName:              v.tName,
			qName:              v.qName,
		}
		samConsensus(s)

		if !fileio.AreEqual(v.outFile, v.outFile_expected) {
			t.Errorf("Error in samConsensus: generating output fa file")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
		if v.multiFaDir != "" {
			for i := range v.multiFaOutAndExpectedFiles {
				if !fileio.AreEqual(i, v.multiFaOutAndExpectedFiles[i]) {
					t.Errorf("Error in samConsensus: output multiFa file did not match expected.")
				} else {
					err = os.Remove(i)
				}
			}
		}
	}
}
