package main

import(
	"log"
	"fmt"
	"flag"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
	"sort"
	)


func assemblyStats(infile string, outfile string) {
	records := fasta.Read(infile)
	var scaf bool = false
	var scafLen int = 0
	var scafList []int
	var genomeLength, N50, sum int

	for i := 0; i < len(records); i++ {
		for k := 0; k < len(records[i].Seq); k++ {
			if scaf {
				if records[i].Seq[k] == dna.N || dna.IsLower(records[i].Seq[k]) {
					scaf = false
					scafList = append(scafList, scafLen)
					scafLen = 0
				} else {
					scafLen++
				}
			} else {
				if !(records[i].Seq[k] == dna.N || dna.IsLower(records[i].Seq[k])) {
					scaf = true
					scafLen++
				}
			}
		}
	}

	for i := 0; i < len(scafList); i++ {
		genomeLength += scafList[i]
	}
	sort.Ints(scafList)
	halfGenome := genomeLength / 2

	for i := len(scafList) -1; i > -1; i-- {
		sum += scafList[i]
		if sum > halfGenome {
			N50 = scafList[i]
			break
		}
	}

	file := fileio.EasyCreate(outfile)
	defer file.Close()

	var err error
	_, err = fmt.Fprintf(file, "Assembly Name: %s\n", infile)
	_, err = fmt.Fprintf(file, "halfGenome: %d\n", halfGenome)
	_, err = fmt.Fprintf(file, "genomeLength: %d\n", genomeLength)
	_, err = fmt.Fprintf(file, "Number of contigs: %d\n", len(scafList))
	_, err = fmt.Fprintf(file, "Largest Contig: %d\n", scafList[len(scafList)-1])
	_, err = fmt.Fprintf(file, "N50: %d\n", N50)
	common.ExitIfError(err)
}

func usage() {
	fmt.Print(
		"assemblyStats: Provides information about the number of scaffolds, including the N50, number of scaffolds, and distribution of lengths of assembled scaffolds.\n" +
		"Usage:\n" +
		"assemblyStats input.fa output.txt\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	assemblyStats(infile, outfile)
}
