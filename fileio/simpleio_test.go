package fileio

import (
	"testing"
)

var testFile string = "testdata/big.fa.gz"

func TestSimpleReader(t *testing.T) {
	answer := ezReaderTest(testFile)
	reader := NewSimpleReader(testFile)
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String() != answer[i] {
			t.Errorf("Error: line did not match easy reader...\n")
		}
		i++
	}
}

func BenchmarkSimpleReader(b *testing.B) {
	reader := NewSimpleReader(testFile)
	b.ReportAllocs()
	b.ResetTimer()
	for _, done := ReadLine(reader); !done; _, done = ReadLine(reader) {
		//Nothing to assign, testing pure reading of the file
	}
}

func ezReaderTest(filename string) []string {
	var answer []string
	er := EasyOpen(filename)
	for line, done := EasyNextLine(er); !done; line, done = EasyNextLine(er) {
		answer = append(answer, line)
	}
	return answer
}
