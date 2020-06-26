package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"

	"strings"
)

func usage() {
	fmt.Print(
		"faTools - quickly manipulate fasta sequences\n" +
			"Usage:\n" +
			" faTools [options] input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var sortMethod *string = flag.String("sort", "", "sort fasta sequences either `byname or bylen`")
	var faTrim *bool = flag.Bool("trim", false, "return a fragment trimmed by providing `start end` coordinates")
	var start *int = flag.Int("start", 0, "starting base postion for fasta trim zero based, left closed")
	var end *int = flag.Int("end", 0, "ending base postion for fasta trim right open")
	var rename *string = flag.String("rename", "", "provide a prefix to edit the names fasta records")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	inputFile, outputFile := flag.Arg(0), flag.Arg(1)
	switch true {
	case *faTrim:
		//TODO: figure out what is better, create new tags for start and end or just have them as extra arguments
		//decided to go with args because it will thow an error if user does not input an int whereas I might have to catch if it is passed in as a string
		faSubsetFrag(inputFile, outputFile, *start, *end)
	case *sortMethod != "":
		faSortRecords(inputFile, outputFile, *sortMethod)
	case *rename != "":
		setNewPrefix(inputFile, outputFile, *rename)
	default:
		flag.Usage()
	}
}

func faSortRecords(inFile string, outFile string, method string) {
	fa := fasta.Read(inFile)

	file := fileio.EasyCreate(outFile)
	defer file.Close()

	if strings.Contains(method, "byname") {
		fasta.SortByName(fa)
	} else if strings.Contains(method, "bylen") {
		fasta.SortBySeqLen(fa)
	} else {
		log.Fatalf("Error: Method to sort fasta sequences not recognized. Try byname or bylen")
	}
	fasta.WriteToFileHandle(file.File, fa, 50)
}

func faSubsetFrag(inputFile string, outputFile string, start int, end int) {
	faSeq := fasta.Read(inputFile)
	writer := fileio.EasyCreate(outputFile)
	if len(faSeq) != 1 {
		log.Fatalf("Error: Must provide a single fasta sequence for fragment trimming")
	}
	fasta.WriteToFileHandle(writer, []*fasta.Fasta{TrimFasta(faSeq[0], start, end)}, 50)
}

func setNewPrefix(inputFile string, outputFile string, prefix string) {
	faPipe := fasta.NewPipe()
	go faPipe.ReadToChan(inputFile)

	var index int = 0
	writer := fileio.EasyCreate(outputFile)
	defer writer.Close()

	for eachFa := range faPipe.Stream {
		fasta.ChangePrefix(eachFa, prefix, index)
		fasta.WriteFasta(writer, eachFa, 50)
		index++
	}
}

//Trim fasta records by giving start and end coords
func TrimFasta(fa *fasta.Fasta, start int, end int) *fasta.Fasta {
	fa.Seq = fa.Seq[start:end]
	return fa
}
