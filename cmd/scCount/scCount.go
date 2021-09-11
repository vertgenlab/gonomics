// Command Group: "SAM Tools"

//TODO: Input-normalization.

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"log"
	"sort"
	"strings"
)

func scCount(s Settings) {
	var sc sam.SingleCellAlignment
	var err error
	readChan, _ := sam.GoReadToChan(s.InFile)
	genes := gtf.Read(s.GeneFile)
	var geneIntervals []interval.Interval
	var geneIdSlice []string = make([]string, 0)
	var geneIndex = make(map[string]int)
	var headerString string = "Bx"

	for g := range genes {
		geneIdSlice = append(geneIdSlice, g)
	}
	sort.Strings(geneIdSlice)
	for c, g := range geneIdSlice { //c for count, g for gene
		headerString += fmt.Sprintf("\t%s", g)
		geneIndex[genes[g].GeneID] = c
		geneIntervals = append(geneIntervals, genes[g])
		c++
	}

	out := fileio.EasyCreate(s.OutFile)
	_, err = fmt.Fprintf(out, "%s\n", headerString)
	exception.PanicOnErr(err)
	tree := interval.BuildTree(geneIntervals)
	var overlap []interval.Interval
	var currRow Row
	var currGene string
	var firstTime bool = true
	var normFlag bool = false

	var normalizationMap map[string]float64 = make(map[string]float64, 0)
	if s.ExpNormalizationFile != "" {
		normFlag = true
		var line string
		var words []string
		var doneReading bool = false
		exp := fileio.EasyOpen(s.ExpNormalizationFile)
		for line, doneReading = fileio.EasyNextRealLine(exp); !doneReading; line, doneReading = fileio.EasyNextRealLine(exp) {
			words = strings.Split(line, "\t")
			if len(words) != 2 {
				log.Fatalf("Expression normalization input file must be a tab-separated file with two columns per line.")
			}
			normalizationMap[words[0]] = common.StringToFloat64(words[1])
		}
		err = exp.Close()
		exception.PanicOnErr(err)
	}

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
				printRow(out, currRow, normalizationMap, geneIdSlice, normFlag)
			}
			currRow = Row{dna.BasesToString(sc.Bx), make([]float64, len(geneIndex))}
		}
		currRow.Counts[geneIndex[currGene]]++ //increment count for gene in currLine
	}
	printRow(out, currRow, normalizationMap, geneIdSlice, normFlag) //print the last cell
	err = out.Close()
	exception.PanicOnErr(err)
}

func getGeneName(g *gtf.Gene) string {
	return g.GeneID
}

func printRow(out io.Writer, r Row, normalizationMap map[string]float64, geneIdSlice []string, normFlag bool) {
	if normFlag { //if the user passed in a normalization file, we normalize the count for each gene.
		var ok bool
		var val float64
		for i := range r.Counts {
			if val, ok = normalizationMap[geneIdSlice[i]]; ok { //if the current gene is in the normalization map.
				r.Counts[i] = r.Counts[i] * val
			}
		}
	}
	var countString string = fmt.Sprintf("%g", r.Counts[0])
	for i := 1; i < len(r.Counts); i++ {
		countString = countString + fmt.Sprintf("\t%g", r.Counts[i])
	}
	_, err := fmt.Fprintf(out, "%s\t%s\n", r.Bx, countString)
	exception.PanicOnErr(err)
}

//Each Row represents one line of the output tsv, which includes the count for each gene from a particular cell. As counts can be weighted by input normalization factors, counts are represented as floats.
type Row struct {
	Bx     string
	Counts []float64
}

func usage() {
	fmt.Print(
		"scCount - Generate count matrix from single-cell sequencing data.\n\n" +
			"Accepts sam reads aligned to a reference genome that have first been processed with fastqFormat -singleCell -collapseUmi and sorted with mergeSort -singleCellBx\n" +
			"Usage:\n" +
			"scCount reads.sam genes.gtf out.csv\n" +
			"scCount also accepts a tab-separated optional input to declare expression normalization multipliers for each gene.\n" +
			"An expression normalization file must have two columns: the first containing geneIDs, and the second containing a normalization multiplier that can be parsed as a float.\n" +
			"An example expression normalization file is shown below:\n" +
			"\tgene1\t1.2\n" +
			"\tgene2\t1.1\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile               string
	GeneFile             string
	OutFile              string
	ExpNormalizationFile string
}

func main() {
	var expNormFile *string = flag.String("expNormalizationFile", "", "Filename for input-normalization.")
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	s := Settings{
		InFile:               flag.Arg(0),
		GeneFile:             flag.Arg(1),
		OutFile:              flag.Arg(2),
		ExpNormalizationFile: *expNormFile,
	}
	scCount(s)
}
