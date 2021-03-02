package main

import (
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fasta"
)

func usage() {
	fmt.Print(
		"")
	flag.PrintDefaults()
}
func main() {
	var expectedNumArgs int = 2	
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	file := flag.Arg(0)
	fq := fastq.GoReadToChan(file)
	out := fileio.EasyCreate(flag.Arg(1))
	defer out.Close()
	for i := range fq {
		fa := &fasta.Fasta{Name: i.Name, Seq: i.Seq}
		fasta.WriteFasta(out, fa, 50)
	}
}
