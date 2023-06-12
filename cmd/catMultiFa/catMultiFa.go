// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
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
			log.Fatalf("Each file to be concatenated must have the same number of entries. Expected %d, found %d in the file named %s.\n", len(ans), len(curr), fileList[i])
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
			"OR\n" +
			"catMultiFa -list file.list" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var minNumArgs int = 2
	var o *string = flag.String("o", "stdout", "Redirect the output to a user-given filename.")
	var list *string = flag.String("list", "", "User may supply a line-delimited list of files to concatenate instead of providing them as arguments.")
	var lineLength *int = flag.Int("lineLength", 50, "Specify the line length in the output file.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *list != "" {
		if len(flag.Args()) > 0 {
			log.Fatalf("catMultiFa accepts either files as arguments or a list of files with the -list option. However, it does not support the simultaneous usage of arguments and a file list.")
		}
	} else if len(flag.Args()) < minNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting at least %d arguments, but got %d\n", minNumArgs, len(flag.Args()))
	}
	var fileList []string
	if *list == "" {
		fileList = flag.Args()
	} else {
		fileList = fileio.Read(*list)
	}
	catMultiFa(fileList, *o, *lineLength)
}
