package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

var BedSplitTests = []struct {
	Mode        string
	InFile      string
	OutDir      string
	GzipOut     bool
	ExpectedDir string
}{
	{Mode: "byName",
		InFile:      "testdata/test.bed",
		OutDir:      "testdata/tmpByName",
		GzipOut:     false,
		ExpectedDir: "testdata/expectedByName"},
	{Mode: "byName",
		InFile:      "testdata/test.bed",
		OutDir:      "testdata/tmpGzipByName",
		GzipOut:     true,
		ExpectedDir: "testdata/expectedGzipByName"},
	{Mode: "byChrom",
		InFile:      "testdata/test.bed",
		OutDir:      "testdata/tmpByChrom",
		GzipOut:     false,
		ExpectedDir: "testdata/expectedByChrom"},
	{Mode: "byChrom",
		InFile:      "testdata/test.bed",
		OutDir:      "testdata/tmpGzipByChrom",
		GzipOut:     true,
		ExpectedDir: "testdata/expectedGzipByChrom"},
}

func TestBedSplit(t *testing.T) {
	var s Settings
	var outFiles, expectedFiles, currOutFileContents, currExpectedFileContents []string
	var currOutFile string
	var currFileIndex, currLineIndex int
	var err error
	var testsFailed bool = false
	for _, v := range BedSplitTests {
		testsFailed = false
		s = Settings{
			Mode:    v.Mode,
			InFile:  v.InFile,
			OutDir:  v.OutDir,
			GzipOut: v.GzipOut,
		}
		bedSplit(s)
		outFiles, err = filepath.Glob(filepath.Join(v.OutDir, "*"))
		exception.PanicOnErr(err)
		expectedFiles, err = filepath.Glob(filepath.Join(v.ExpectedDir, "*"))
		exception.PanicOnErr(err)
		if len(outFiles) != len(expectedFiles) {
			t.Errorf("Error: in bedSplit, number of output files was not as expected.\n")
			testsFailed = true
		}

		if !testsFailed {
			for currFileIndex, currOutFile = range outFiles {
				currOutFileContents = fileio.Read(currOutFile)
				currExpectedFileContents = fileio.Read(expectedFiles[currFileIndex])
				if strings.Compare(filepath.Base(currOutFile), filepath.Base(expectedFiles[currFileIndex])) != 0 {
					t.Errorf("Error: in bedSplit, output file name was not as expected.\n")
					testsFailed = true
				}
				if len(currOutFileContents) != len(currExpectedFileContents) {
					t.Errorf("Error: in bedSplit, number of lines in output file: %s was not as expected.\n", currOutFile)
					testsFailed = true
				}
				for currLineIndex = range currOutFileContents {
					if strings.Compare(currOutFileContents[currLineIndex], currExpectedFileContents[currLineIndex]) != 0 {
						t.Errorf("Error: in bedsplit, output file did not have expected contents. Offending line: %s\n", currOutFileContents[currLineIndex])
						testsFailed = true
					}
				}
			}
		}

		if !testsFailed {
			err = os.RemoveAll(v.OutDir)
			exception.PanicOnErr(err)
		}
	}
}
