package fileio

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"strings"
	"testing"

	"github.com/klauspost/pgzip"
	"github.com/vertgenlab/gonomics/exception"
)

// example going straight to file.
func writeDnaFile(file *os.File) {
	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	var err error
	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
	exception.PanicOnErr(err)
}

// example using io.Writer.
func writeDnaIo(file io.Writer) {
	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	var err error
	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
	exception.PanicOnErr(err)
}

// example using fileio.EasyWriter.
func writeDnaFileio(file *EasyWriter) {
	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	var err error
	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
	exception.PanicOnErr(err)
}

// used to uncompress the file in the repo.
func copyFile(inputFilename string, outputFilename string) {
	var er *EasyReader
	var ew *EasyWriter
	var line string
	var done bool
	var err error

	er = EasyOpen(inputFilename)
	defer er.Close()
	ew = EasyCreate(outputFilename)
	defer ew.Close()

	for line, done = EasyNextLine(er); !done; line, done = EasyNextLine(er) {
		_, err = fmt.Fprintf(ew, "%s\n", line)
		exception.PanicOnErr(err)
	}
}

func BenchmarkRead(b *testing.B) {
	unzip := "testdata/big.fa"
	copyFile("testdata/big.fa.gz", unzip)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		var er *EasyReader
		var done bool

		er = EasyOpen(unzip)

		for _, done = EasyNextLine(er); !done; _, done = EasyNextLine(er) {
		}
		er.Close()
	}
	EasyRemove(unzip)
}

func BenchmarkReadGz(b *testing.B) {
	for n := 0; n < b.N; n++ {
		var er *EasyReader
		var done bool

		er = EasyOpen("testdata/big.fa.gz")

		for _, done = EasyNextLine(er); !done; _, done = EasyNextLine(er) {
		}
		er.Close()
	}
}

func BenchmarkWriteFileio(b *testing.B) {
	test := "testdata/testWrite.dna"
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate(test)

		for i := 0; i < 10000; i++ {
			writeDnaFileio(ew)
		}
		ew.Close()
	}
	EasyRemove(test)
}

func BenchmarkWriteFileioGz(b *testing.B) {
	test := "testdata/testWrite.dna.gz"
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate(test)

		for i := 0; i < 10000; i++ {
			writeDnaFileio(ew)
		}
		ew.Close()
	}
	EasyRemove(test)
}

func BenchmarkWriteIo(b *testing.B) {
	test := "testdata/testWrite.dna"
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate(test)

		for i := 0; i < 10000; i++ {
			writeDnaIo(ew)
		}
		ew.Close()
	}
	EasyRemove(test)
}

func BenchmarkWriteIoGz(b *testing.B) {
	test := "testdata/testWrite.dna.gz"
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate(test)

		for i := 0; i < 10000; i++ {
			writeDnaIo(ew)
		}
		ew.Close()
	}
	EasyRemove(test)
}

func BenchmarkWriteFile(b *testing.B) {
	test := "testdata/testWrite.dna"
	for n := 0; n < b.N; n++ {
		osf, err := os.Create(test)
		exception.PanicOnErr(err)

		for i := 0; i < 10000; i++ {
			writeDnaFile(osf)
		}
		osf.Close()
	}
	EasyRemove(test)
}

var testfile2 string = "testdata/smallTest"
var line4 string = "#shhhh this line is a secret"
var line5 string = "Hello World"
var line6 string = "I am a gopher"

func TestRead(t *testing.T) {
	file := Read(testfile2)
	if file[0] != line5 {
		t.Errorf("problem with read")
	}
	if file[1] != line6 {
		t.Errorf("problem with read")
	}
}

func TestWrite(t *testing.T) {
	var createdFile string
	var err error
	createdFile = "test.txt"
	var fileContent []string = []string{line4, line5, line6}
	Write(createdFile, fileContent)
	if !AreEqual(testfile, createdFile) {
		t.Errorf("problem with fileio.Write()")
	} else {
		err = os.Remove("test.txt")
		exception.PanicOnErr(err)
	}
}

func TestStdinGzip(t *testing.T) {
	testCases := []string{
		"fileio stdin with gzip",
		"stdin.gz",
	}

	for _, input := range testCases {
		// Create a temporary file to simulate stdin
		mockStdin, err := os.CreateTemp("", "mockStdin.gz")
		if err != nil {
			t.Fatalf("Error: Failed to create temporary file: %v", err)
		}
		defer mockStdin.Close()
		defer EasyRemove(mockStdin.Name())

		// Gzip the input data
		var gzippedInput bytes.Buffer
		gzWriter := pgzip.NewWriter(&gzippedInput)
		if _, err := gzWriter.Write([]byte(input)); err != nil {
			t.Fatalf("Error: Failed to gzip input: %v", err)
		}
		if err := gzWriter.Close(); err != nil {
			t.Fatalf("Error: Failed to close gzip writer: %v", err)
		}

		// Write the gzipped input to the mock stdin
		if _, err := io.Copy(mockStdin, &gzippedInput); err != nil {
			t.Fatalf("Error: Failed to write simulated input: %v", err)
		}

		// Seek to the beginning of the mock stdin
		if _, err := mockStdin.Seek(0, io.SeekStart); err != nil {
			t.Fatalf("Error: Failed to seek to the start: %v", err)
		}

		// Set os.Stdin to the mock stdin
		os.Stdin = mockStdin

		// Call the function being tested
		file := EasyOpen("/dev/stdin")
		defer file.Close()

		// Read from the file and compare with the original input
		var reader strings.Builder
		if _, err := io.Copy(&reader, file); err != nil {
			t.Fatalf("Error: Failed to read from simulated stdin: %v", err)
		}
		if reader.String() != input {
			t.Errorf("Error: Mismatch for input %q. Expected: %q, got %q", input, input, reader.String())
		}
	}
}
