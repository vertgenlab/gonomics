package fileio

import (
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"sync"
	"testing"
)

var bufsizes = []int{
	0, 7, 16, 23, 32, 46, 64, 93, 128, 1024,
}

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
	unzip := "testdata/big.fa"
	copyFile("testdata/big.fa.gz", unzip)
	b.ReportAllocs()
	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		reader := NewByteReader(unzip)
		for _, done := ReadLine(reader); !done; _, done = ReadLine(reader) {
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
		for _, done := ReadLine(reader); !done; _, done = ReadLine(reader) {
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

func TestWriter(t *testing.T) {
	var data [8192]byte

	for i := 0; i < len(data); i++ {
		data[i] = byte(' ' + i%('~'-' '))
	}
	w := new(bytes.Buffer)
	for i := 0; i < len(bufsizes); i++ {
		for j := 0; j < len(bufsizes); j++ {
			nwrite := bufsizes[i]
			bs := bufsizes[j]

			// Write nwrite bytes using buffer size bs.
			// Check that the right amount makes it out
			// and that the data is correct.

			w.Reset()
			buf := NewByteWriterSize(w, bs)
			context := fmt.Sprintf("nwrite=%d bufsize=%d", nwrite, bs)
			n, e1 := buf.Write(data[0:nwrite])
			if e1 != nil || n != nwrite {
				t.Errorf("%s: buf.Write = %d, %v", context, n, e1)
				continue
			}
			if e := buf.Flush(); e != nil {
				t.Errorf("%s: buf.Flush = %v", context, e)
			}

			written := w.Bytes()
			if len(written) != nwrite {
				t.Errorf("%s: %d bytes written", context, len(written))
			}
			for l := 0; l < len(written); l++ {
				if written[l] != data[l] {
					t.Errorf("wrong bytes written")
					t.Errorf("want=%q", data[0:len(written)])
					t.Errorf("have=%q", written)
				}
			}
		}
	}
}

func TestEmptyWrite(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, defaultBufSize)

	n, err := writer.Write([]byte{})
	if err != nil {
		t.Fatalf("Write failed: %v", err)
	}
	if n != 0 {
		t.Fatalf("Expected 0 bytes written, got %d", n)
	}
	writer.Flush()
}

func TestFlushWithoutWrite(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, defaultBufSize)

	err := writer.Flush()
	if err != nil {
		t.Fatalf("Flush failed: %v", err)
	}
}

func TestLargeWrite(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, 10) // Small buffer to force flush

	data := []byte("This is a very large write that exceeds the buffer size.")
	n, err := writer.Write(data)
	if err != nil {
		t.Fatalf("Write failed: %v", err)
	}
	if n != len(data) {
		t.Fatalf("Expected %d bytes written, got %d", len(data), n)
	}
	writer.Flush()

	if buf.String() != string(data) {
		t.Fatalf("Expected %s, got %s", string(data), buf.String())
	}
}

func TestConcurrentWrites(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, 10) // Small buffer to force flush

	data := []byte("Concurrent write.\n")
	var wg sync.WaitGroup

	for i := 0; i < 10; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			_, err := writer.Write(data)
			if err != nil {
				t.Errorf("Write failed: %v", err)
			}
		}()
	}

	wg.Wait()
	writer.Flush()

	expected := bytes.Repeat(data, 10)
	if buf.String() != string(expected) {
		t.Fatalf("Expected %s, got %s", string(expected), buf.String())
	}
}

func TestResetWriter(t *testing.T) {
	var buf1, buf2 bytes.Buffer
	writer := NewByteWriterSize(&buf1, defaultBufSize)

	data := []byte("Initial write.")
	writer.Write(data)
	writer.Flush()

	writer.Reset(&buf2)
	data = []byte("After reset.")
	writer.Write(data)
	writer.Flush()

	if buf1.String() != "Initial write." {
		t.Fatalf("Expected buf1 to have 'Initial write.', got %s", buf1.String())
	}
	if buf2.String() != "After reset." {
		t.Fatalf("Expected buf2 to have 'After reset.', got %s", buf2.String())
	}
}

func TestCloseWithoutWrite(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, defaultBufSize)

	err := writer.Close()
	if err != nil {
		t.Fatalf("Close failed: %v", err)
	}
}

func TestWriteAfterClose(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, defaultBufSize)
	writer.Close()

	_, err := writer.Write([]byte("Write after close."))
	if !errors.Is(err, io.ErrClosedPipe) {
		t.Fatalf("Expected io.ErrClosedPipe, got %v", err)
	}
}
func TestMultipleFlushCalls(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, defaultBufSize)

	data := []byte("Multiple flushes.")
	writer.Write(data)
	writer.Flush()
	writer.Flush()
	writer.Flush()

	if buf.String() != string(data) {
		t.Fatalf("Expected %s, got %s", string(data), buf.String())
	}
}

func TestInterleavedWritesAndFlushes(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, defaultBufSize)

	data1 := []byte("First write.")
	data2 := []byte("Second write.")
	writer.Write(data1)
	writer.Flush()
	writer.Write(data2)
	writer.Flush()

	expected := "First write.Second write."
	if buf.String() != expected {
		t.Fatalf("Expected %s, got %s", expected, buf.String())
	}
}

func TestWriteStringWithLongString(t *testing.T) {
	var buf bytes.Buffer
	writer := NewByteWriterSize(&buf, 10) // Small buffer to force flush

	data := "This is a very large string that exceeds the buffer size."
	n, err := writer.WriteString(data)
	if err != nil {
		t.Fatalf("WriteString failed: %v", err)
	}
	if n != len(data) {
		t.Fatalf("Expected %d bytes written, got %d", len(data), n)
	}
	writer.Flush()

	if buf.String() != data {
		t.Fatalf("Expected %s, got %s", data, buf.String())
	}
}

// TestGzipWrite tests writing data to a gzip file using ByteWriter.
func TestGzipWrite(t *testing.T) {
	filename := "test_output.gz"
	writer := NewByteWriter(filename)

	// Data to be written
	data := "This is a test line.\n"

	// Write data to the writer
	_, err := writer.Write([]byte(data))
	if err != nil {
		t.Fatalf("Write failed: %v", err)
	}
	// Flush the writer to ensure all data is written
	err = writer.Flush()
	if err != nil {
		t.Fatalf("Flush failed: %v", err)
	}

	// Close the writer
	err = writer.Close()
	if err != nil {
		t.Fatalf("Close failed: %v", err)
	}

	// Read the gzip file and verify its contents
	file, err := os.Open(filename)
	if err != nil {
		t.Fatalf("Failed to open gzip file: %v", err)
	}
	defer file.Close()
	defer os.Remove(filename)

	gz, err := gzip.NewReader(file)
	if err != nil {
		t.Fatalf("Failed to create gzip reader: %v", err)
	}
	defer gz.Close()

	content, err := ioutil.ReadAll(gz)
	if err != nil {
		t.Fatalf("Failed to read gzip content: %v", err)
	}

	if string(content) != data {
		t.Fatalf("Expected %q, got %q", data, string(content))
	}
}
func BenchmarkEasyWriter(b *testing.B) {
	filename := "benchmark_easywriter_test.txt"
	data := []byte("This is a benchmark test line.\n") // Data to write
	writer := EasyCreate(filename)

	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time

	for i := 0; i < b.N; i++ {
		for j := 0; j < 12; j++ {
			writer.Write(data)
			writer.Close()
		}
	}
	// Cleanup the benchmark file
	os.Remove(filename)
}

func BenchmarkBytekWriter(b *testing.B) {
	filename := "benchmark_fasterwriter_test.txt"
	file, err := os.Create(filename)
	if err != nil {
		b.Fatalf("Failed to create file: %v", err)
	}
	defer os.Remove(filename)

	data := []byte("This is a benchmark test line.\n")
	bw := NewByteWriterSize(file, defaultBufSize)

	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time
	for i := 0; i < b.N; i++ {
		// Write multiple lines within the benchmark iteration
		for j := 0; j < 12; j++ {
			_, err := bw.Write(data)
			if err != nil {
				b.Fatalf("Write failed: %v", err)
			}
		}
	}
	bw.Flush() // Ensure final flush for accurate measurement
}
