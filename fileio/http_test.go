package fileio

import (
	"testing"
)

var httpTestingLookUp = []struct {
	local    string
	internet string
}{
	{"../README.md", "https://raw.githubusercontent.com/vertgenlab/gonomics/main/README.md"},
}

// TestHttpReader will parse a local text file and a file uploaded to an internew server and compare the lines
// an check if they are exact matches.
func TestHttpReader(t *testing.T) {
	for _, test := range httpTestingLookUp {
		// ReadUrl is the function we are testing to see the reader is processing the lines of data correctly
		internetData := ReadUrl(test.internet)

		// reading local file
		localReader := NewSimpleReader(test.local)
		defer localReader.Close()
		var index int = 0
		// loop to read inputed local text file
		for line, done := ReadLine(localReader); !done; line, done = ReadLine(localReader) {
			if line.String() == internetData[index] {
				t.Logf("%s\n", internetData[index])
				index++
			}
		}
		if index != len(internetData) {
			t.Errorf("Error: http did not read the same number of lines as the local file...\n")
		}
	}
}
