package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
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
	var sortMethod *string = flag.String("sort", "", "sort fasta sequences either `[byname/bylen]`")
	var faTrim *string = flag.String("trim", "", "return a fragment trimmed by providing `start:end` coordinates, zero based, left closed, right open")
	var rename *string = flag.String("rename", "", "provide a `prefix` to edit the names fasta records")

	//need to keep out file separate from flags for merging fasta for *.fa in a directory
	var merge *bool = flag.Bool("merge", false, "takes a list of .fasta files to merge into one file `*.fa or [file.fa file.fa ...]`")
	var outMerge *string = flag.String("out", "/dev/stdout", "output filename for merging `.fa`sta records")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	inputFile, outputFile := flag.Arg(0), flag.Arg(1)
	switch true {
	case *faTrim != "":
		faSubsetFrag(inputFile, outputFile, *faTrim)
	case *sortMethod != "":
		faSortRecords(inputFile, outputFile, *sortMethod)
	case *rename != "":
		setNewPrefix(inputFile, outputFile, *rename)
	case *merge:
		combineFastaFiles(flag.Args(), *outMerge)
	default:
		flag.Usage()
	}
}

func combineFastaFiles(files []string, outputFile string) {
	ans := fasta.NewPipe()
	ans.Wg.Add(1)
	go fasta.ReadMultiFilesToChan(ans, files)
	go fasta.WritingChannel(outputFile, ans.StdOut, ans.Wg)
	ans.Wg.Wait()
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

func faSubsetFrag(inputFile string, outputFile string, window string) {
	faSeq := fasta.Read(inputFile)
	writer := fileio.EasyCreate(outputFile)
	pos := strings.Split(window, ":")
	if len(faSeq) != 1 {
		log.Fatalf("Error: Must provide a single fasta sequence for fragment trimming")
	}
	fasta.WriteToFileHandle(writer, []*fasta.Fasta{TrimFasta(faSeq[0], common.StringToInt(pos[0]), common.StringToInt(pos[1]))}, 50)
}

func setNewPrefix(inputFile string, outputFile string, prefix string) {
	faPipe := fasta.NewPipe()
	go fasta.ReadToChan(inputFile, faPipe.StdOut)

	var index int = 0
	writer := fileio.EasyCreate(outputFile)
	defer writer.Close()

	for eachFa := range faPipe.StdOut {
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
