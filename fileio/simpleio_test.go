package fileio

import (
	"testing"
)

func TestSimpleReader(t *testing.T) {
	answer := ezReaderTest("testdata/big.fa.gz")
	reader := NewSimpleReader("testdata/big.fa.gz")
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String() != answer[i] {
			t.Errorf("Error: line did not match easy reader...\n")
		}
		i++
	}
}

func TestLineExceedsDefaultBufferSize(t *testing.T) {
	answer := ezReaderTest("testdata/longlined.vcf")
	reader := NewSimpleReader("testdata/longlined.vcf")
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String() != answer[i] {
			t.Errorf("Error: line did not match easy reader...\n")
		}
		i++
	}
}

func BenchmarkSimpleReader(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewSimpleReader("testdata/big.fa")
		for _, done := ReadLine(reader); !done; _, done = ReadLine(reader) {
			//Nothing to assign, testing pure reading of the file
		}
	}
}

func BenchmarkSimpleReaderGz(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewSimpleReader("testdata/big.fa.gz")
		for _, done := ReadLine(reader); !done; _, done = ReadLine(reader) {
			//Nothing to assign, testing pure reading of the file
		}
	}
}

func BenchmarkEasyReaderReg(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		var reader *EasyReader
		var done bool
		reader = EasyOpen("testdata/big.fa")
		for _, done = EasyNextLine(reader); !done; _, done = EasyNextLine(reader) {
		}
	}
}

func BenchmarkEasyReaderGz(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		var reader *EasyReader
		var done bool
		reader = EasyOpen("testdata/big.fa.gz")
		for _, done = EasyNextLine(reader); !done; _, done = EasyNextLine(reader) {
		}
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
