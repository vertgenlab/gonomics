package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"log"
	"os"
)

// FilterSettings defines the usage settings for the gtfTools filter subcommand.
type FilterSettings struct {
	InFile           string
	OutFile          string
	GeneNameList     string
	ChromFilter      string
	CodingTranscript bool
}

func filterUsage(filterFlags *flag.FlagSet) {
	fmt.Printf("gtfTools filter - a tool for filtering GTF records.\n" +
		"Usage:\n" +
		"gtfTools filter in.gtf out.gtf\n" +
		"options:\n")
	filterFlags.PrintDefaults()
}

// parseFilterArgs is the main function for the gtfTools filter subcommand. It parses options
// and executes the function gtfFilter.
func parseFilterArgs() {
	var expectedNumArgs int = 2
	var err error
	filterFlags := flag.NewFlagSet("filter", flag.ExitOnError)
	var geneNameList *string = filterFlags.String("geneNameList", "", "Specify a new-line delimited file containing the geneNames for records to be retained.\n")
	var chromFilter *string = filterFlags.String("chromFilter", "", "Specify a chromosome for which all transcript records will be in the output. All transcripts must be on the filtering chromosome in order to pass the filter. Can be used in combination with Gene Name Filter option, all records for a gene must be on the filtering chromosome.\n")
	var codingTranscript *bool = filterFlags.Bool("codingTranscript", false, "In -codingTranscript mode, only transcripts that are coding, i.e. have 1 or more CDS, are retained.\n")
	err = filterFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	filterFlags.Usage = func() { filterUsage(filterFlags) }

	if len(filterFlags.Args()) != expectedNumArgs {
		filterFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(filterFlags.Args()))
	}

	inFile := filterFlags.Arg(0)
	outFile := filterFlags.Arg(1)

	s := FilterSettings{
		InFile:           inFile,
		OutFile:          outFile,
		GeneNameList:     *geneNameList,
		ChromFilter:      *chromFilter,
		CodingTranscript: *codingTranscript,
	}
	gtfFilter(s)
}

func gtfFilter(s FilterSettings) {
	var pass, foundInMap bool
	var err error
	var geneNameMap = make(map[string]bool)
	var currTranscripts []*gtf.Transcript
	out := fileio.EasyCreate(s.OutFile)

	if s.GeneNameList != "" {
		geneNames := fileio.Read(s.GeneNameList)
		for i := range geneNames {
			geneNameMap[geneNames[i]] = true
		}
	}

	records := gtf.Read(s.InFile)
	for currGene := range records { // for each gene
		pass = true

		if s.CodingTranscript { // option to detect coding transcript
			for _, currTranscript := range records[currGene].Transcripts { // for each transcript
				if len(currTranscript.Exons) == 0 { // if transcript has no exon, then transcript is non-coding
					pass = false
				} else { // if transcript has exon(s)
					// if we find at least 1 exon that has 1 CDS, then transcript is coding
					foundCDS := false
					for _, currExon := range currTranscript.Exons { // for each exon
						if currExon.Cds != nil { // if exon has CDS, then found "at least 1 exon that has 1 CDS", transcript is coding
							foundCDS = true
							break
						}
					}
					if !foundCDS { // after checking each exon, if no CDS, then transcript is non-coding
						pass = false
					}
				}
			}
		}

		if s.GeneNameList != "" && s.ChromFilter == "" {
			if _, foundInMap = geneNameMap[records[currGene].GeneName]; !foundInMap {
				pass = false
			}
		} else if s.ChromFilter != "" && s.GeneNameList == "" {
			currTranscripts = records[currGene].Transcripts
			for t := range currTranscripts {
				if currTranscripts[t].Chr != s.ChromFilter {
					pass = false
					break
				}
			}
		} else if s.GeneNameList != "" && s.ChromFilter != "" {
			if _, foundInMap = geneNameMap[records[currGene].GeneName]; !foundInMap {
				currTranscripts = records[currGene].Transcripts
				for t := range currTranscripts {
					if currTranscripts[t].Chr != s.ChromFilter {
						pass = false
						break
					}
				}
			}
		}
		if pass {
			gtf.WriteToFileHandle(out, records[currGene])
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
