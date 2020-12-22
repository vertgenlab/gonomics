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
		internetData := EasyOpen(test.internet)
		var httpData []string
		for data, done := EasyNextLine(internetData); !done; data, done = EasyNextLine(internetData) {
			httpData = append(httpData, data)
		}

		localData := EasyOpen(test.local)
		var index int = 0
		for line, done := EasyNextLine(localData); !done; line, done = EasyNextLine(localData) {
			if line == httpData[index] {
				index++
			} else {
				t.Errorf("Error: fetching data over http did not equal local copy...\n")
			}
		}

	}
}
