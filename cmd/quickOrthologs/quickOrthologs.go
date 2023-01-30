package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
	"log"
)

type Settings struct {
	TName               string
	QName               string
	GtfFile             string
	ChainFile           string
	ChromSizes          string
	OutFile             string
	Unmapped            string
	CanonicalTranscript bool
}

func quickOrthologs(s Settings) {
	var err error
	var currLift lift.Lift
	var overlap []interval.Interval
	var destinationChrom string
	var destinationStart, destinationEnd int
	var chainIntervals []interval.Interval
	var tssBeds []bed.Bed
	genes := gtf.Read(s.GtfFile)
	chroms := chromInfo.ReadToMap(s.ChromSizes)
	if s.CanonicalTranscript {
		tssBeds = gtf.GenesToCanonicalTranscriptsTssBed(genes, chroms)
	} else {
		tssBeds = gtf.GenesToTssBed(genes, chroms)
	}
	chainChan, _ := chain.GoReadToChan(s.ChainFile)
	for val := range chainChan {
		chainIntervals = append(chainIntervals, val)
	}
	tree := interval.BuildTree(chainIntervals)
	out := fileio.EasyCreate(s.OutFile)
	un := fileio.EasyCreate(s.Unmapped)

	_, err = fmt.Fprintf(out, "#geneName\ttName\ttChrom\ttStart\ttEnd\tqName\tqChrom\tqStart\tqEnd\n")

	for _, currTss := range tssBeds {
		if currTss.Name != "" { //we'll discard entries with no geneName for now
			currLift = currTss
			overlap = interval.Query(tree, currLift, "any")
			if len(overlap) > 1 {
				_, err = fmt.Fprintf(un, "Record below maps to multiple chains:\n")
				exception.PanicOnErr(err)
				currTss.WriteToFileHandle(un)
			} else if len(overlap) == 0 {
				_, err = fmt.Fprintf(un, "Record below has no ortholog in new assembly:\n")
				exception.PanicOnErr(err)
				currTss.WriteToFileHandle(un)
			} else {
				destinationChrom, destinationStart, destinationEnd = lift.LiftCoordinatesWithChain(overlap[0].(chain.Chain), currLift)
				_, err = fmt.Fprintf(out, "%s\t%s\t%v\t%v\t%v\t%s\t%s\t%v\t%v\n", currTss.Name, s.TName, currLift.GetChrom(), currLift.GetChromStart(), currLift.GetChromEnd(), s.QName, destinationChrom, destinationStart, destinationEnd)
				exception.PanicOnErr(err)
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
	err = un.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"quickOrthologs - Find corresponding orthologous transcription start" +
			"sites.\n" +
			"Usage:\n" +
			"quickOrthologs tName qName genes.gtf liftOver.chain inputRef.chrom.sizes ortholog.txt unmapped.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var canonicalTranscript *bool = flag.Bool("canonicalTranscript", false, "Use only the canonical transcript for each gene.")
	var expectedNumArgs int = 7

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	tName := flag.Arg(0)
	qName := flag.Arg(1)
	gtfFile := flag.Arg(2)
	chainFile := flag.Arg(3)
	chromSizes := flag.Arg(4)
	outFile := flag.Arg(5)
	unmapped := flag.Arg(6)
	s := Settings{
		TName:               tName,
		QName:               qName,
		GtfFile:             gtfFile,
		ChainFile:           chainFile,
		ChromSizes:          chromSizes,
		OutFile:             outFile,
		Unmapped:            unmapped,
		CanonicalTranscript: *canonicalTranscript,
	}

	quickOrthologs(s)
}
