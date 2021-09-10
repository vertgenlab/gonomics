// Command Group: "SAM Tools"

//TODO: Input-normalization.

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"log"
)

func scCount(s Settings) {
	var sc sam.SingleCellAlignment
	var err error
	readChan, _ := sam.GoReadToChan(s.InFile)
	genes := gtf.Read(s.GeneFile)
	var geneIntervals []interval.Interval
	var geneIndex = make(map[string]int)
	var count int = 0
	var headerString string = "Bx"
	for g := range genes {
		headerString += fmt.Sprintf("\t%s", g)
		geneIndex[genes[g].GeneID] = count
		geneIntervals = append(geneIntervals, genes[g])
		count++
	}
	out := fileio.EasyCreate(s.OutFile)
	_, err = fmt.Fprintf(out,"%s\n", headerString)
	exception.PanicOnErr(err)
	tree := interval.BuildTree(geneIntervals)
	var overlap []interval.Interval
	var currRow Row
	var currGene string
	var firstTime bool = true

	for i := range readChan {
		sc = sam.ToSingleCellAlignment(i)
		overlap = interval.Query(tree, &i, "any")
		if len(overlap) == 0 {
			continue
		}
		if len(overlap) > 1 {
			//TODO: Some transcripts overlap, but we can figure out from exon lines which gene to assign a read to.
			log.Fatalf("The following input SAM record maps to multiple genes in the gtf file:\n%v.", sam.ToString(i))
		}
		currGene = getGeneName(overlap[0].(*gtf.Gene))
		if dna.BasesToString(sc.Bx) != currRow.Bx {
			if firstTime {
				firstTime = false
			} else {
				printRow(out, currRow)
			}
			currRow = Row{dna.BasesToString(sc.Bx), make([]float64, len(geneIndex))}
		}
		currRow.Counts[geneIndex[currGene]]++//increment count for gene in currLine
	}
	printRow(out, currRow)//print the last cell
	err = out.Close()
	exception.PanicOnErr(err)
}

func getGeneName(g *gtf.Gene) string {
	return g.GeneID
}

func printRow(out io.Writer, r Row) {
	var countString string = fmt.Sprintf("%g", r.Counts[0])
	for i := 1; i < len(r.Counts); i++ {
		countString = countString + fmt.Sprintf("\t%g", r.Counts[i])
	}
	_, err := fmt.Fprintf(out,"%s\t%s\n", r.Bx, countString)
	exception.PanicOnErr(err)
}

//Each Row represents one line of the output tsv, which includes the count for each gene from a particular cell. As counts can be weighted by input normalization factors, counts are represented as floats.
type Row struct {
	Bx	string
	Counts []float64
}

func usage() {
	fmt.Print(
		"scCount - Generate count matrix from single-cell sequencing data." +
			"Accepts sam reads aligned to a reference genome that have first been processed with fastqFormat -singleCell -collapseUmi and sorted with mergeSort -singleCellBx\n" +
			"Usage:\n" +
			"scCount reads.sam genes.gtf out.csv\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile string
	GeneFile string
	OutFile string
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		InFile: flag.Arg(0),
		GeneFile: flag.Arg(1),
		OutFile: flag.Arg(2),
	}

	scCount(s)
}