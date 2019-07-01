package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"log"
)

func bedAnnotate(in string, rec string, outfile string, field *int, prefix string) {
	var records []*bed.Bed = bed.Read(rec)
	var infile []*bed.Bed = bed.Read(in)

	for i := 0; i < len(infile); i++ {
		for j := 0; j < len(records); j++ {
			if overlap(infile[i], records[j]) && !contains(infile[i].Annotation, prefix + records[j].Name) {
				infile[i].Annotation = append(infile[i].Annotation, prefix + records[j].Name)
			}
		}
	}

	bed.Write(outfile, infile, 7)
}

func usage() {
	fmt.Print(
		"bedAnnotate - Appends a field from bed2 into overlapping segments of bed1\n" +
		"Usage:\n" +
		"bedAnnotate bed1.bed records.bed output.bed\n" +
		"options:\n")
	flag.PrintDefaults()
}

func overlap(bed1 *bed.Bed, bed2 *bed.Bed) bool {
	if bed1.Chrom != bed2.Chrom {
		return false
	} else if bed1.ChromEnd < bed2.ChromStart || bed2.ChromEnd < bed1.ChromStart {
		return false
	}
	return true
}

func contains(s []string, e string) bool {
    for _, a := range s {
        if a == e {
            return true
        }
    }
    return false
}

func main() {
	var expectedNumArgs int = 3
	var field *int = flag.Int("field", 6, "Specifies the bed field to add annotation")
	var prefix *string = flag.String("prefix", "", "Adds a prefix to the annotation value")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	records := flag.Arg(1)
	outfile := flag.Arg(2)

	bedAnnotate(infile, records, outfile, field, *prefix)
}
