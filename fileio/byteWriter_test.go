package fileio

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"sync"
	"testing"
)

var bufsizes = []int{
	0, 7, 16, 23, 32, 46, 64, 93, 128, 1024,
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
			buf := NewWriterSize(w, bs)
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

func BenchmarkEasyWriter(b *testing.B) {
	filename := "benchmark_easywriter_test.txt"
	data := []byte("This is a benchmark test line.\n") // Data to write
	writer := EasyCreate(filename)

	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time

	for i := 0; i < b.N; i++ {
		for j := 0; j < 100; j++ {
			writer.Write(data)
			writer.Close()
		}
	}
	// Cleanup the benchmark file
	os.Remove(filename)
}

func BenchmarkFasterkWriter(b *testing.B) {
	filename := "benchmark_fasterwriter_test.txt"
	file := MustCreate(filename)

	data := []byte("This is a benchmark test line.\n")
	bw := NewWriter(file)
	writer := NewWriter(bw)

	b.ReportAllocs() // Enable memory allocation reporting
	b.ResetTimer()   // Reset timer to exclude setup time
	for i := 0; i < b.N; i++ {
		// Write multiple lines within the benchmark iteration
		for j := 0; j < 100; j++ {
			writer.Write(data)
		}
		writer.Flush()
		// Flush less often for a fairer comparison (tune this)
	}
	writer.Flush() // Ensure final flush for accurate measurement
	os.Remove(filename)
}

func BenchmarkLargeConcurrentWrite(b *testing.B) {
	bw := NewWriterSize(io.Discard, 1024) // Use io.Discard to avoid file I/O overhead

	data := []byte("This is a benchmark test line.\n")

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var wg sync.WaitGroup

		// Write multiple lines within the benchmark iteration
		for j := 0; j < 100; j++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				_, err := bw.Write(data)
				if err != nil {
					b.Fatal(err)
				}
			}()
		}

		// Wait for all writes in this iteration to complete
		wg.Wait()

		// Flush the buffer only once per benchmark iteration
		err := bw.Flush()
		if err != nil {
			b.Fatal(err)
		}
	}
}

// func TestCorrectNessByteWriter(t *testing.T) {
// 	file := MustCreate("testdata/bytewriter.test")
// 	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

// 	defer file.Close()
// 	var err error
// 	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
// 	exception.PanicOnErr(err)

// 	buf := new(bytes.Buffer)
//     bw := NewWriter(buf)
//     bw.WriteString(fakeDna)

// 	if buf.String() != fakeDna {
//         t.Errorf("Expected %q, got %q", fakeDna, buf.String())
//     }
// }
