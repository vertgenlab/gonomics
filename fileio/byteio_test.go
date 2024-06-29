package fileio

import (
	"fmt"
	"os"
	"sync"
	"testing"
)

var bufsizes = []int{
	0, 7, 16, 23, 32, 46, 64, 93, 128, 1024,
}

func TestByteioReader(t *testing.T) {
	answer := ezReaderTest("testdata/big.fa")
	reader := NewByteReader("testdata/big.fa")
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String() != answer[i] {
			t.Errorf("Error: line did not match easy reader...\n")
		} else {
			fmt.Printf("%q == %s\n", line, answer[i])
			i++
		}
	}
}

func TestLineExceedsDefaultBufferSize(t *testing.T) {
	answer := ezReaderTest("testdata/lineExceedsBufferSize.txt")
	reader := NewByteReader("testdata/lineExceedsBufferSize.txt")
	var i int = 0
	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
		if line.String()== answer[i] {
			fmt.Printf("%q == %s\n", line, answer[i])
			i++
		} else {
			t.Errorf("Error: byteio reader did not process bytes into buffer beyond the default buffer size %q != %s ...\n", line, answer[i])
		}
	}
}

func BenchmarkByteioReaderReg(b *testing.B) {
	unzip := "testdata/big.fa"
	copyFile("testdata/big.fa.gz", unzip)
	b.ReportAllocs()
	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		reader := NewByteReader(unzip)
		for _, done :=ReadLine(reader); !done; _, done = ReadLine(reader) {
			//Nothing to assign, testing pure reading of the file
		}
	}
	EasyRemove(unzip)
}

func BenchmarkByteioReaderGz(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewByteReader("testdata/big.fa.gz")
		for _, done := NextByteLine(reader); !done; _, done = NextByteLine(reader) {
			//Nothing to assign, testing pure reading of the file
		}
	}
}

func BenchmarkEasyReaderReg(b *testing.B) {
	unzip := "testdata/big.fa"
	copyFile("testdata/big.fa.gz", unzip)
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		var reader *EasyReader
		var done bool
		reader = EasyOpen(unzip)
		for _, done = EasyNextLine(reader); !done; _, done = EasyNextLine(reader) {
		}
	}
	EasyRemove(unzip)
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

// func TestByteWriterGz(t *testing.T) {
// 	reader := NewByteReader("testdata/big.fa.gz")
// 	gz := "testdata/byte.gz"
// 	defer reader.Close()

// 	// Read and store data from the original GZ file
// 	var originalData [][]byte
// 	for line, done := ReadLine(reader); !done; line, done = ReadLine(reader) {
// 		originalData = append(originalData, line.Bytes())
// 	}

// 	// Write the data to a new GZ file
// 	writer := NewByteWriter(gz)
// 	for _, line := range originalData {
// 		_, err := writer.Write(line)
// 		exception.PanicOnErr(err)
// 	}
// 	writer.Close() // Ensure the GZ file is properly closed

// 	// Read the data back from the newly created GZ file
// 	result := NewByteReader(gz)
// 	defer result.Close()
// 	var writtenData [][]byte

// 	for line, done := ReadLine(result); !done; line, done = ReadLine(result) {
// 		writtenData = append(writtenData, line.Bytes())
// 	}

// 	if len(originalData) == len(writtenData) {
// 		for i := 0; i < len(writtenData); i++ {
// 			if string(originalData[i]) != string(writtenData[i]) {
// 				t.Logf("Error: %q != %q", originalData[i], writtenData[i])
// 			}
// 		}
// 	} else {
// 		t.Logf("Error: size of results do not match %d != %d", len(originalData), len(writtenData))
// 	}
// }

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


// func TestWriteStringWithLongString(t *testing.T) {
// 	var buf bytes.Buffer
// 	writer := NewByteWriterSize(&buf, 10) // Small buffer to force flush

// 	data := "This is a very large string that exceeds the buffer size."
// 	n, err := writer.WriteString(data)
// 	if err != nil {
// 		t.Fatalf("WriteString failed: %v", err)
// 	}
// 	if n != len(data) {
// 		t.Fatalf("Expected %d bytes written, got %d", len(data), n)
// 	}
// 	writer.Flush()

// 	if buf.String() != data {
// 		t.Fatalf("Expected %s, got %s", data, buf.String())
// 	}
// }

// TestGzipWrite tests writing data to a gzip file using ByteWriter.
// func TestGzipWrite(t *testing.T) {
// 	filename := "test_output.gz"
// 	writer := NewByteWriter(filename)

// 	// Data to be written
// 	data := "This is a test line.\n"

// 	// Write data to the writer
// 	_, err := writer.Write([]byte(data))
// 	if err != nil {
// 		t.Fatalf("Write failed: %v", err)
// 	}
// 	// Flush the writer to ensure all data is written
// 	err = writer.Flush()
// 	if err != nil {
// 		t.Fatalf("Flush failed: %v", err)
// 	}

// 	// Close the writer
// 	err = writer.Close()
// 	if err != nil {
// 		t.Fatalf("Close failed: %v", err)
// 	}

// 	// Read the gzip file and verify its contents
// 	file, err := os.Open(filename)
// 	if err != nil {
// 		t.Fatalf("Failed to open gzip file: %v", err)
// 	}
// 	defer file.Close()
// 	defer os.Remove(filename)

// 	gz, err := gzip.NewReader(file)
// 	if err != nil {
// 		t.Fatalf("Failed to create gzip reader: %v", err)
// 	}
// 	defer gz.Close()

// 	content, err := io.ReadAll(gz)
// 	if err != nil {
// 		t.Fatalf("Failed to read gzip content: %v", err)
// 	}

// 	if string(content) != data {
// 		t.Fatalf("Expected %q, got %q", data, string(content))
// 	}
// }
func BenchmarkEasyWriterReg(b *testing.B) {
	filename := "benchmark_easywriter_test.txt"
	data := []byte("This is a benchmark test line.\n") // Data to write
	//var err error
	writer := EasyCreate(filename)
	defer os.Remove(filename)

	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time

	for i := 0; i < b.N; i++ {
		writer.Write(data)
	}
}

func BenchmarkByteWriterReg(b *testing.B) {
    filename := "benchmark_fasterwriter_test.txt"
    writer := NewByteWriter(filename)

    defer os.Remove(filename)

    data := []byte("This is a benchmark test line.\n")

    b.ReportAllocs() // Enable memory allocation reporting
    b.ResetTimer()   // Reset timer to exclude setup time

    for i := 0; i < b.N; i++ {
		writer.Write(data)
    }

}

func BenchmarkEasyWriterGz(b *testing.B) {
	filename := "benchmark_easywriter_test.txt.gz"
	data := []byte("This is a benchmark test line.\n") // Data to write
	//var err error
	writer := EasyCreate(filename)
	

	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time

	for i := 0; i < b.N; i++ {
		writer.Write(data)
	}
	os.Remove(filename)
}
func BenchmarkBytekWriterGz(b *testing.B) {
	filename := "benchmark_fasterwriter_test.txt.gz"
	bw := NewByteWriter(filename)
	var wg sync.WaitGroup
	

	data := []byte("This is a benchmark test line.\n")
	var j int
	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time
	for i := 0; i < b.N/4; i++ {
		wg.Add(4)
		for j =0; j < 4;j++ {
			go func(index int) {
				bw.Write(data)
				if index%4==0 {
					bw.Flush()
				}
				wg.Done()
			}(j)
		}
		wg.Wait()
	}
    os.Remove(filename)
}

// func BenchmarkSmallConcurrentWrite(b *testing.B) {
// 	b.ReportAllocs()
// 	bw := NewByteWriterSize(io.Discard, 1024)
// 	data := []byte("This is a benchmark test line.\n")

// 	var wg sync.WaitGroup
// 	for i := 0; i < b.N; i++ {
// 		wg.Add(12)
// 		for j :=0; j < 12;j++ {
// 			go func(index int) {
// 				bw.Write(data)
// 				if index%33 == 0 {
// 					bw.Flush()
// 				}
// 				wg.Done()
// 			}(j)
// 		}
// 		wg.Wait()
// 	}
// }