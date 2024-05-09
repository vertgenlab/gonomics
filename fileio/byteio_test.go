package fileio

import (
	"testing"
)

func TestByteioReader(t *testing.T) {
	answer := ezReaderTest("testdata/big.fa.gz")
	reader := NewByteReader("testdata/big.fa.gz")
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String() != answer[i] {
			t.Errorf("Error: line did not match easy reader...\n")
		}
		i++
	}
}

func TestLineExceedsDefaultBufferSize(t *testing.T) {
	answer := ezReaderTest("testdata/lineExceedsBufferSize.txt")
	reader := NewByteReader("testdata/lineExceedsBufferSize.txt")
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String() == answer[i] {
			i++
		} else {
			t.Errorf("Error: byteio reader did not process bytes into buffer beyond the default buffer size %s != %s ...\n", line.String(), answer[i])
		}
	}
}

func BenchmarkByteioReader(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewByteReader("testdata/big.fa")
		for _, done := ReadLine(reader); !done; _, done = ReadLine(reader) {
			//Nothing to assign, testing pure reading of the file
		}
	}
}

func BenchmarkByteioReaderGz(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewByteReader("testdata/big.fa.gz")
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

func TestStringToIntSlice(t *testing.T) {
	data := "0,1,2,3,4,5,6,7,8,9,"
	expected := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	answer := StringToIntSlice(data)
	if len(answer) != 10 {
		t.Errorf("Error: problem converting string to int slice")
	}
	for i := range answer {
		if answer[i] != expected[i] {
			t.Errorf("Error: problem converting string to int slice, %d != %d...\n", answer[i], expected[i])
		}
	}
}

func TestIntSliceToString(t *testing.T) {
	expected := "0,1,2,3,4,5,6,7,8,9,"
	data := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	answer := IntSliceToString(data)
	if answer != expected {
		t.Errorf("Error: problem converting int slice to string")
	}
}
