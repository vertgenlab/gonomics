// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"unicode"
)

func usage() {
	fmt.Print(
		"faMaskIupac - N-mask extended IUPAC nucleotides (includes UWSMKRYBDHV)\n\n" +
			"Usage:\n" +
			"  faMaskIupac [options] in.fa\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func main() {
	var output *string = flag.String("o", "stdout", "Output fasta")
	flag.Parse()
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	if len(flag.Args()) != 1 {
		usage()
		return
	}

	infile := flag.Arg(0)

	faMaskIupac(infile, *output)
}

func faMaskIupac(input, output string) {
	var err error
	out := fileio.EasyCreate(output)
	in := fileio.EasyOpen(input)

	var line string
	var doneReading bool
	for line, doneReading = fileio.EasyNextLine(in); !doneReading; line, doneReading = fileio.EasyNextLine(in) {
		if strings.HasPrefix(line, "#") || strings.HasPrefix(line, ">") {
			fileio.WriteToFileHandle(out, line)
		}
		fileio.WriteToFileHandle(out, mask([]byte(line)))
	}

	err = in.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func mask(line []byte) string {
	for i := range line {
		switch line[i] {
		case 'U', 'u':
			fallthrough
		case 'W', 'w':
			fallthrough
		case 'S', 's':
			fallthrough
		case 'M', 'm':
			fallthrough
		case 'K', 'k':
			fallthrough
		case 'R', 'r':
			fallthrough
		case 'Y', 'y':
			fallthrough
		case 'B', 'b':
			fallthrough
		case 'D', 'd':
			fallthrough
		case 'H', 'h':
			fallthrough
		case 'V', 'v':
			if unicode.IsUpper(rune(line[i])) {
				line[i] = 'N'
			} else {
				line[i] = 'n'
			}
		}
	}
	return string(line)
}
