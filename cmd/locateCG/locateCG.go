package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

type Settings struct {
	InFa      string //input fasta file with 2 aligned sequences to compare
	ChromName string //e.g. chr1
	Outfile   string //output file name
	Compare   bool   //if inputting an aligned fasta file with 2 genomes to compare
	CGtype    string //"gain", "loss", "cons", "all"
}

type CGdiff struct {
	Chrom    string //derived from ChromName
	StartPos int    //start pos of CG site in reference genome of sequence 1
	EndPos   int    //end pos of CG site in reference genome of sequence 1
	Type     string //gain or loss in sequence 1 relative to sequence 2
	Ref      string //sequence 1 base(s)
	Alt      string //sequence 2 base(s)
	AlnStart int    //corresponding alignment position of reference StartPos
	AlnEnd   int    //corresponding alignment position of reference EndPos
}

// helper function for file writing in --compare mode
func toFile(fileName string, records []CGdiff) {
	var err error

	//create file
	file := fileio.EasyCreate(fileName)

	//write header
	fileio.WriteToFileHandle(file, "Chrom\tRefStart\tRefEnd\tType\tRef\tAlt\tAlnStart\tAlnEnd")

	//loop through structs and write each struct to a line
	for _, r := range records {
		//1 added to end coordinates when writing to file for consistency with BED formatting
		line := fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d", r.Chrom, r.StartPos, r.EndPos+1, r.Type, r.Ref, r.Alt, r.AlnStart, r.AlnEnd+1)
		fileio.WriteToFileHandle(file, line)
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

// -compare mode check if constant CG
func isCons(firstSeq1, firstSeq2, secondSeq1, secondSeq2 dna.Base) bool {
	if firstSeq1 == dna.C && firstSeq2 == dna.G && secondSeq1 == dna.C && secondSeq2 == dna.G {
		return true
	} else {
		return false
	}
}

// -compare mode check if gained CG in seq1
func isGain(firstSeq1, firstSeq2, secondSeq1, secondSeq2 dna.Base) bool {
	if firstSeq1 == dna.C && firstSeq2 == dna.G && !(secondSeq1 == dna.C && secondSeq2 == dna.G) {
		return true
	} else {
		return false
	}
}

// -compare mode check if lost CG in seq1
func isLoss(firstSeq1, firstSeq2, secondSeq1, secondSeq2 dna.Base) bool {
	if !(firstSeq1 == dna.C && firstSeq2 == dna.G) && secondSeq1 == dna.C && secondSeq2 == dna.G {
		return true
	} else {
		return false
	}
}

// default function for one fasta sequence
func locateCG(s Settings) {
	var base1, base2 dna.Base
	var output []bed.Bed
	f := fasta.Read(s.InFa)

	//check that there are not more than 1 fasta sequence
	if len(f) == 2 {
		log.Fatalf("Error: expecting exactly one record in fasta file, but got 2. Are you trying to compare two sequences? If yes, use --compare mode.\n")
	}

	if len(f) > 2 {
		log.Fatalf("Error: expecting exactly one record in fasta file, but got %d. Using --compare mode will allow for comparison of 2 sequences, but not more than 2.\n", len(f))
	}

	seq := f[0].Seq
	//if sequence is empty
	if len(seq) == 0 {
		log.Fatalf("Error: fasta sequence is empty.\n")
	}

	//for length of fasta sequence
	for i := 0; i < len(seq)-1; i++ {
		//store dna.Bases in pairs
		base1, base2 = seq[i], seq[i+1]

		//if "CG" present
		if base1 == dna.C && base2 == dna.G {
			output = append(output, bed.Bed{Chrom: s.ChromName, ChromStart: i, ChromEnd: i + 2, FieldsInitialized: 3})
		}
	}

	//create output file name and initialize file, write header line
	outfileName := s.Outfile
	bed.Write(outfileName, output)
	fmt.Println("CG sites found and written to", outfileName)

}

// --compare mode function
func compareCG(s Settings) {
	var f1, f2, s1, s2 dna.Base
	var refBases, altBases string
	var output []CGdiff

	//read in fasta file
	f := fasta.Read(s.InFa)

	//check that there are exactly 2 sequences to compare
	if len(f) != 2 {
		log.Fatalf("Error: --compare mode expects exactly two sequences in fasta, but got %d\n", len(f))
	}

	//store sequence 1 and 2
	firstSeq := f[0].Seq
	secondSeq := f[1].Seq

	//if any sequences empty
	if len(firstSeq) == 0 || len(secondSeq) == 0 {
		log.Fatalf("Missing or empty sequences for seq1 or seq2.")
	}

	//check for equal sequence lengths
	if len(firstSeq) != len(secondSeq) {
		log.Fatalf("Seq1 and seq2 not equal in length.")
	}

	//pre-allocate size of output
	//output = make([]CGsite, 0, 10000)

	//initialize pair of matched reference/alignment positions
	refStart := 0
	alnStart := 0

	//for length of first sequence
	for i := 0; i < len(firstSeq)-1; i++ {

		//store dna.Bases in pairs in both sequences
		f1, f2 = firstSeq[i], firstSeq[i+1]
		s1, s2 = secondSeq[i], secondSeq[i+1]

		//skip if not bases (gaps)
		if !(dna.DefineBase(f1)) || !(dna.DefineBase(f2)) || !(dna.DefineBase(s1)) || !(dna.DefineBase(s2)) {
			continue
		}

		//store "Ref" bases from sequence 1 and "Alt" bases from sequence 2
		refBases = dna.BaseToString(f1) + dna.BaseToString(f2)
		altBases = dna.BaseToString(s1) + dna.BaseToString(s2)

		switch s.CGtype {
		case "cons":
			if isCons(f1, f2, s1, s2) {
				startPos := fasta.AlnPosToRefPosCounter(f[0], i, refStart, alnStart)
				endPos := startPos + 1
				output = append(output, CGdiff{Chrom: s.ChromName, StartPos: startPos, EndPos: endPos, Type: "cons", Ref: refBases, Alt: altBases, AlnStart: i, AlnEnd: i + 1})
			}
		case "gain":
			if isGain(f1, f2, s1, s2) {
				startPos := fasta.AlnPosToRefPosCounter(f[0], i, refStart, alnStart)
				endPos := startPos + 1
				output = append(output, CGdiff{Chrom: s.ChromName, StartPos: startPos, EndPos: endPos, Type: "gain", Ref: refBases, Alt: altBases, AlnStart: i, AlnEnd: i + 1})
			}
		case "loss":
			if isLoss(f1, f2, s1, s2) {
				startPos := fasta.AlnPosToRefPosCounter(f[0], i, refStart, alnStart)
				endPos := startPos + 1
				output = append(output, CGdiff{Chrom: s.ChromName, StartPos: startPos, EndPos: endPos, Type: "loss", Ref: refBases, Alt: altBases, AlnStart: i, AlnEnd: i + 1})
			}
		default:
			log.Fatalf("Unknown CpG comparison type: %s. Options: gain, loss, cons", s.CGtype)
		}
		//if a CG site has already been recorded in "output"
		if len(output) > 0 {
			//reset the refStart to the last reference position
			refStart = output[len(output)-1].StartPos
			//reset the alnStart to the last corresponding alignment position
			alnStart = output[len(output)-1].AlnStart
			//this will now be the new ref/aln position pair for conversion
		}
	}

	//write to file
	outfileName := s.Outfile
	toFile(outfileName, output)
	fmt.Println("CG comparisons found and written to", outfileName)
}

func usage() {
	fmt.Print("locateCG chr*.fa.gz chromName outfile\n" +
		"\tRecords positional information of CpG sites. If inputting a single sequence/genome fasta,\n" +
		"\treports position of all CpG sites. If inputting fasta of 2 aligned genomes, use -compare mode \n" +
		"\tand option -cgtype to identify CpG site changes in sequence 1 relative to sequence 2. \n" +
		"\t All CpG site comparisons are reported from substitutions only, not indels. \n" +
		"\t \n" +
		"\tARGS:\n" +
		"\tchr*fa.gz - input fasta file (per chromosome)\n" +
		"\tchromName - name of chromosome e.g. chr1\n" +
		"\toutfile - name of output file to write to (.bed for single-genome, .txt for -compare mode)\n" +
		"\t \n" +
		"\tUSAGE:\n" +
		"\tlocateCG chr*fa.fa.gz chromName outfile.bed -> one genome\n" +
		"\tOR\n" +
		"\tlocateCG -compare -cgtype={'gain', 'loss', 'cons'} chr*fa.gz chromName outfile.txt -> two genomes\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	//expect 3 args: input fasta file, chromosome name, and outfile name
	var expectedNumArgs int = 3
	var compare *bool = flag.Bool("compare", false, "Will compare 2 aligned genomes from an alignment file. Must use with -cgtype to specify type of CpG comparison to report.")
	var cgtype *string = flag.String("cgtype", "", "Type(s) of CpG site comparisons to report when comparing genomes. Options: 'gain', 'loss', 'cons'. All comparisons are reported for sequence 1 relative to sequence 2.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	//if number of provided args is not 3
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		InFa:      flag.Arg(0),
		ChromName: flag.Arg(1),
		Outfile:   flag.Arg(2),
		Compare:   *compare,
		CGtype:    *cgtype,
	}

	//run locate() function on provided args; will output a slice of structs
	if *compare {
		if *cgtype == "" {
			log.Fatalf("Error: Must specify --cgtype in --compare mode: 'gain', 'loss', 'cons'.")
		}
		if *cgtype != "gain" && *cgtype != "loss" && *cgtype != "cons" {
			log.Fatalf("Error: Unrecognized CG type. Options: 'gain', 'loss', 'cons'.")
		}
		compareCG(s)
	} else {
		locateCG(s)
	}
}
