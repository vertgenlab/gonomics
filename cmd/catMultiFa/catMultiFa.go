// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func catMultiFa(fileList []string, o string, lineLength int) {
	if len(fileList) < 1 {
		log.Fatalf("Must provide at least one file to cat. File list is empty.")
	}
	var err error
	ans := fasta.Read(fileList[0])
	var curr []fasta.Fasta

	for i := 1; i < len(fileList); i++ {
		curr = fasta.Read(fileList[i])
		if len(curr) != len(ans) {
			log.Fatalf("Each file to be concatenated must have the same number of entries. Expected %v, found %v in the file named %s.\n", len(ans), len(curr), fileList[i])
		}
		for j := range curr {
			if curr[j].Name != ans[j].Name {
				log.Fatalf("Each file to be concatenated must contain the same names in the same order. Expected %v, found %v in the file named %s.\n", ans[j].Name, curr[j].Name, fileList[i])
			}
			ans[j].Seq = append(ans[j].Seq, curr[j].Seq...)
		}
	}

	out := fileio.EasyCreate(o)
	fasta.WriteToFileHandle(out, ans, lineLength)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"catMultiFa - Concatenate multiFa alignment files by sequence. Accepts an arbitrary number of input fasta files. " +
			"Input multiFa files must have the same number and order of records. Additional sequence from each file is appended onto the sequence of each entry." +
			"For large number of files, use the list option to avoid the command character limit. Default to standard output.\n" +
			"Usage:\n" +
			"catMultiFa multi.fa ...\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var minNumArgs int = 2
	var o *string = flag.String("o", "stdout", "Redirect the output to a user-given filename.")
	var list *string = flag.String("list", "", "")
	var lineLength *int = flag.Int("lineLength", 50, "Specify the linelength in the output file.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *list != "" {
		minNumArgs = 0
	}

	if len(flag.Args()) < minNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting at least %d arguments, but got %d\n", minNumArgs, len(flag.Args()))
	}
	var fileList []string
	var err error
	if *list == "" {
		fileList = make([]string, len(flag.Args()))
		for i := range flag.Args() {
			fileList[i] = flag.Arg(i)
		}
	} else {
		file := fileio.EasyOpen(*list)
		fileList = make([]string, 0)
		var line string
		var doneReading bool
		for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
			fileList = append(fileList, line)
		}
		err = file.Close()
		exception.PanicOnErr(err)
	}
	catMultiFa(fileList, *o, *lineLength)
}
