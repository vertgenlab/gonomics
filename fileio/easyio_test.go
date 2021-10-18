package fileio

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"os"
	"testing"
)

// example going straight to file
func writeDnaFile(file *os.File) {
	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	var err error
	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
	exception.PanicOnErr(err)
}

// example using io.Writer
func writeDnaIo(file io.Writer) {
	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	var err error
	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
	exception.PanicOnErr(err)
}

// example using fileio.EasyWriter
func writeDnaFileio(file *EasyWriter) {
	var fakeDna string = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	var err error
	_, err = fmt.Fprintf(file, "%s\n", fakeDna)
	exception.PanicOnErr(err)
}

// used to uncompress the file in the repo
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
	copyFile("testdata/big.fa.gz", "testdata/big.fa")
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		var er *EasyReader
		var done bool

		er = EasyOpen("testdata/big.fa")

		for _, done = EasyNextLine(er); !done; _, done = EasyNextLine(er) {
		}
		er.Close()
	}
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
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate("testdata/testWrite.dna")

		for i := 0; i < 10000; i++ {
			writeDnaFileio(ew)
		}
		ew.Close()
	}
}

func BenchmarkWriteFileioGz(b *testing.B) {
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate("testdata/testWrite.dna.gz")

		for i := 0; i < 10000; i++ {
			writeDnaFileio(ew)
		}
		ew.Close()
	}
}

func BenchmarkWriteIo(b *testing.B) {
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate("testdata/testWrite.dna")

		for i := 0; i < 10000; i++ {
			writeDnaIo(ew)
		}
		ew.Close()
	}
}

func BenchmarkWriteIoGz(b *testing.B) {
	for n := 0; n < b.N; n++ {
		var ew *EasyWriter = EasyCreate("testdata/testWrite.dna.gz")

		for i := 0; i < 10000; i++ {
			writeDnaIo(ew)
		}
		ew.Close()
	}
}

func BenchmarkWriteFile(b *testing.B) {
	for n := 0; n < b.N; n++ {
		osf, err := os.Create("testdata/testWrite.dna")
		exception.PanicOnErr(err)

		for i := 0; i < 10000; i++ {
			writeDnaFile(osf)
		}
		osf.Close()
	}
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

