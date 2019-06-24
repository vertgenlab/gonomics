package wig

import (
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/in_test.wig"},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		_ = Read(test.filename)
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Wig
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)

		if !AllEqual(Read(tempFile), Read("testdata/in_test.wig")) {
			t.Errorf("Read and write Wig files were not the same")
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}

/*
func test_wig(inFile string, outFile string) {
	var wigList []wig.Wig
	wigList, err = wig.Read(inFile)
	if err != nil {
		return err
	}

	wig.Write(outFile, wigList)
}

func usage() {
	fmt.Print(
		"wig_test - Tests wig package functions\n" +
			"Usage:\n" +
			" wig_test input.wig output.wig\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	test_wig(inFile, outFile)
}
*/
