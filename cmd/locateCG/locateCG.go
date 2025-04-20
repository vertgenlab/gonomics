package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

type Settings struct {
	InFa      string //input fasta file with 2 aligned sequences to compare
	ChromName string //e.g. chr1
}

type CGsite struct {
	Chrom    string //derived from ChromName
	StartPos int    //start pos of CG site in reference genome of sequence 1
	EndPos   int    //end pos of CG site in reference genome of sequence 1
	Type     string //gain or loss in sequence 1 relative to sequence 2
	AlnStart int    //corresponding alignment position of reference StartPos
	AlnEnd   int    //corresponding alignment position of reference EndPos
	Ref      string //sequence 1 base(s)
	Alt      string //sequence 2 base(s)
}

func locateCG(s Settings) {
	var err error
	var output []CGsite
	var f1, f2, s1, s2 dna.Base
	var refBases, altBases string

	//read in fasta file
	f := fasta.Read(s.InFa)

	if len(f) > 2 {
		log.Fatalf("Error: expecting exactly two records in multiFa, but got %d\n", len(f))
	}

	//specify sequence 1 and 2
	firstSeq := f[0].Seq
	secondSeq := f[1].Seq

	//if no sequences found
	if len(firstSeq) == 0 || len(secondSeq) == 0 {
		log.Fatalf("Missing or empty sequences for seq1 or seq2.")
	}
	//add check for equal lengths
	if len(firstSeq) != len(secondSeq) {
		log.Fatalf("Seq1 and seq2 not equal in length.")
	}
	//pre-allocate size of output
	//output = make([]CGsite, 0, 10000)

	//initialize pair of matched reference/alignment positions
	refStart := 0
	alnStart := 0

	//fmt.Printf("Beginning loop\n")
	//for length of first sequence
	for i := 0; i < len(firstSeq)-1; i++ {
		//store dna.Bases
		f1, f2 = firstSeq[i], firstSeq[i+1]
		s1, s2 = secondSeq[i], secondSeq[i+1]
		if !(dna.DefineBase(f1)) || !(dna.DefineBase(f2)) || !(dna.DefineBase(s1)) || !(dna.DefineBase(s2)) {
			continue
		}
		refBases = ""
		altBases = ""
		if f1 != s1 {
			refBases += dna.BaseToString(f1)
			altBases += dna.BaseToString(s1)
		}
		if f2 != s2 {
			refBases += dna.BaseToString(f2)
			altBases += dna.BaseToString(s2)
		}

		//if there is a CG in sequence 1
		if f1 == dna.C && f2 == dna.G {
			//but no CG in sequence 2
			if !(s1 == dna.C && s2 == dna.G) {
				//convert alignment position to reference position
				startPos := fasta.AlnPosToRefPosCounter(f[0], i, refStart, alnStart)
				endPos := startPos + 1
				//add to "gained CpGs"
				output = append(output, CGsite{s.ChromName, startPos, endPos, "gain", i, (i + 1), refBases, altBases})
			}
			//if there is no CG in sequence 1 but there is a CG in sequence 2
		} else if s1 == dna.C && s2 == dna.G {
			//convert alignment position to reference position
			startPos := fasta.AlnPosToRefPosCounter(f[0], i, refStart, alnStart)
			endPos := startPos + 1
			//add to "lost CpGs"
			output = append(output, CGsite{s.ChromName, startPos, endPos, "loss", i, (i + 1), refBases, altBases})
		}

		//if a CG gain or loss has already been recorded in "output"
		if len(output) > 0 {
			//reset the refStart to the last reference position
			refStart = output[len(output)-1].StartPos
			//reset the alnStart to the last corresponding alignment position
			alnStart = output[len(output)-1].AlnStart
			//this will now be the new ref/aln position pair for conversion
		}
	}
	//fmt.Printf("Finished loop\n")

	//create output file name and initialize file, write header line
	outfileName := s.ChromName + "CGs.txt"
	outfile := fileio.EasyCreate(outfileName)
	fileio.WriteToFileHandle(outfile, "Chrom\tRefStart\tRefEnd\tType\tAlnStart\tAlnEnd\tRef\tAlt")

	//write each struct to a line of the output file
	for _, p := range output {
		line := fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s", p.Chrom, p.StartPos, p.EndPos, p.Type, p.AlnStart, p.AlnEnd, p.Ref, p.Alt)
		fileio.WriteToFileHandle(outfile, line)
	}

	err = outfile.Close()
	exception.PanicOnErr(err)

	fmt.Println("CG gains & losses found and written to", outfileName)
}

func usage() {
	fmt.Print("locateCG chr*.fa.gz chromName \n" +
		"\tRecords positional information of CpG site changes in sequence 1 relative\n" +
		"\tto sequence 2 in alignment fasta. Only records changes resulting from substitutions, not indels.\n" +
		"\tARGS:\n" +
		"\tchr*fa.gz - input fasta file  with 2 aligned sequences to compare\n" +
		"\tchromName - name of chromosome e.g. chr1\n")
}

func main() {
	//expect 2 args: input fasta file and chromosome name
	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	//if number of provided args is not 2
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		InFa:      flag.Arg(0),
		ChromName: flag.Arg(1),
	}

	//run locate() function on provided args; will output a slice of structs
	locateCG(s)
}
