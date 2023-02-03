// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

type Settings struct {
	SamFileName string
	RefFile string
	OutFile string
	VcfFile string
	ChainFile string
	SubstitutionsOnly bool
	InsertionThreshold float64
}

const bufferSize = 10_000_000

func samConsensus(s Settings) {
	var err error
	var outVcfFile *fileio.EasyWriter
	var currConsensus sam.Consensus
	var currChrom string
	var answerPos int//position in the currChrom of the answer
	var firstTime bool = true
	var emptyRoomInBuffer, refPos, i int
	var newBufferRoom []dna.Base = make([]dna.Base, bufferSize)
	var currFaIndex int
	var positionsToSkip int
	var refAllele []dna.Base

	if s.InsertionThreshold < 0 || s.InsertionThreshold > 1 {
		log.Fatalf("InsertionThreshold option must be a value between 0 and 1. Found: %v.\n", s.InsertionThreshold)
	}

	ref := fasta.Read(s.RefFile)
	for i = range ref {
		dna.AllToLower(ref[i].Seq)
	}
	refMap := fasta.ToMap(ref)

	reads, header := sam.GoReadToChan(s.SamFileName)
	piles := sam.GoPileup(reads, header, false, nil, nil)
	if s.VcfFile != "" {
		outVcfFile = fileio.EasyCreate(s.VcfFile)
		_, err = fmt.Fprintf(outVcfFile, "%s\n", strings.Join(vcf.NewHeader(s.SamFileName).Text, "\n"))
		exception.PanicOnErr(err)
	}

	var answer []fasta.Fasta = make([]fasta.Fasta, len(ref))
	for i = range ref {
		answer[i] = fasta.Fasta{Name: ref[i].Name, Seq: make([]dna.Base, bufferSize)}
	}

	for p := range piles {
		if positionsToSkip > 0 {
			positionsToSkip--
			refPos++
			continue
		}
		if firstTime {
			currChrom = header.Chroms[p.RefIdx].Name
			currFaIndex = getIndexForName(answer, currChrom)
			//fmt.Printf("Length of answer[currChrom]: %v.\n", len(answer[currFaIndex].Seq))
			emptyRoomInBuffer = bufferSize
			answerPos, refPos = 0, 0
			firstTime = false
		}
		if currChrom != header.Chroms[p.RefIdx].Name {//if we've moved onto a new chromosome.
			for refPos < len(refMap[currChrom]) {
				answer[currFaIndex].Seq[answerPos] = refMap[currChrom][refPos]
				refPos++
				answerPos++
				emptyRoomInBuffer--
				if emptyRoomInBuffer < 1 {
					answer[currFaIndex].Seq = append(answer[currFaIndex].Seq, newBufferRoom...)
					emptyRoomInBuffer += bufferSize
				}
			}
			answer[currFaIndex].Seq = answer[currFaIndex].Seq[:len(answer[currFaIndex].Seq) - emptyRoomInBuffer]//clear out empty buffer positions
			//now we set up the new chromosome
			currChrom = header.Chroms[p.RefIdx].Name
			currFaIndex = getIndexForName(answer, currChrom)
			emptyRoomInBuffer = bufferSize
			answerPos, refPos = 0, 0
		}

		//catch up to the current pile position, handles reference positions with no Pile coverage.
		for refPos < int(p.Pos - 1) {
			answer[currFaIndex].Seq[answerPos] = refMap[currChrom][refPos]//refMap is lower case
			emptyRoomInBuffer--
			answerPos++
			refPos++
			if emptyRoomInBuffer < 1 {
				answer[currFaIndex].Seq = append(answer[currFaIndex].Seq, newBufferRoom...)
				emptyRoomInBuffer += bufferSize
			}
		}

		currConsensus = sam.PileConsensus(p, s.SubstitutionsOnly, s.InsertionThreshold)

		//now refPos should equal p.Pos - 1, because of our for loop before
		if refPos != int(p.Pos - 1) {
			log.Fatalf("Something went wrong. RefPos is not equal to p.Pos -1.")
		}
		switch currConsensus.Type {
		case sam.Undefined:
			answer[currFaIndex].Seq[answerPos] = refMap[currChrom][refPos]//refMap is lowercase so we'll get lower case in the answer
			emptyRoomInBuffer--
			if emptyRoomInBuffer < 1 {
				answer[currFaIndex].Seq = append(answer[currFaIndex].Seq, newBufferRoom...)
				emptyRoomInBuffer += bufferSize
			}
			answerPos++
			refPos++
		case sam.Base:
			answer[currFaIndex].Seq[answerPos] = currConsensus.Base
			if s.VcfFile != "" && currConsensus.Base != dna.ToUpper(refMap[currChrom][refPos]) {
				_, err = fmt.Fprintf(outVcfFile, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", currChrom, int64(p.Pos), ".", dna.BaseToString(dna.ToUpper(refMap[currChrom][refPos])), dna.BaseToString(currConsensus.Base), ".", ".", ".", ".")
				exception.PanicOnErr(err)
			}
			emptyRoomInBuffer--
			if emptyRoomInBuffer < 1 {
				answer[currFaIndex].Seq = append(answer[currFaIndex].Seq, newBufferRoom...)
				emptyRoomInBuffer += bufferSize
			}
			answerPos++
			refPos++
		case sam.Insertion:
			answer[currFaIndex].Seq[answerPos] = currConsensus.Base // real reference base at insertion is before the inserted sequence.
			emptyRoomInBuffer--
			if emptyRoomInBuffer < 1 {
				answer[currFaIndex].Seq = append(answer[currFaIndex].Seq, newBufferRoom...)
				emptyRoomInBuffer += bufferSize
			}
			answerPos++
			for i = range currConsensus.Insertion {
				answer[currFaIndex].Seq[answerPos] = currConsensus.Insertion[i]
				emptyRoomInBuffer--
				if emptyRoomInBuffer < 1 {
					answer[currFaIndex].Seq = append(answer[currFaIndex].Seq, newBufferRoom...)
					emptyRoomInBuffer += bufferSize
				}
				answerPos++
			}
			if s.VcfFile != "" {
				_, err = fmt.Fprintf(outVcfFile, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", currChrom, int(p.Pos), ".", dna.BaseToString(dna.ToUpper(refMap[currChrom][refPos])), dna.BasesToString(currConsensus.Insertion), ".", ".", ".", ".")
				exception.PanicOnErr(err)
			}
			refPos++
		case sam.Deletion://this pile corresponds to the deleted position, so nothing is written here, we skip for the duration of the deletion.
			positionsToSkip = currConsensus.Deletion - 1
			if s.VcfFile != "" {
				refAllele = make([]dna.Base, currConsensus.Deletion + 1) // make a dna.Base slice the length of the deletion + 1 (to include base before deletion)
				for i = range refAllele {
					refAllele[i] = dna.ToUpper(refMap[currChrom][refPos + i])//starts with refPos (p.Pos - 1) to include the position before the deletion in the VCF
				}
				_, err = fmt.Fprintf(outVcfFile, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", currChrom, int(p.Pos), ".", dna.BasesToString(refAllele), dna.BaseToString(dna.ToUpper(refMap[currChrom][refPos])), ".", ".", ".", ".")
				exception.PanicOnErr(err)
			}
			refPos++
		}
	}
	// once we're done with the piles we have to add the trailing ref bases and clear the buffer for the last chrom
	for refPos < len(refMap[currChrom]) {
		answer[currFaIndex].Seq[answerPos] = refMap[currChrom][refPos]
		refPos++
		answerPos++
		emptyRoomInBuffer--
	}
	answer[currFaIndex].Seq = answer[currFaIndex].Seq[:len(answer[currFaIndex].Seq) - emptyRoomInBuffer]

	if s.VcfFile != "" {
		err = outVcfFile.Close()
		exception.PanicOnErr(err)
	}

	fasta.Write(s.OutFile, answer)
}

func getIndexForName(f []fasta.Fasta, name string) int {
	for i := range f {
		if f[i].Name == name {
			return i
		}
	}
	log.Fatalf("Name: %s not found in fasta.", name)
	return -1
}

func usage() {
	fmt.Print(
		"samConsensus - Generates a fasta file from a sam over a reference sequence.\n" +
			"Uncovered sequences are converted to lowercase reference sequences.\n" +
			"Usage:\n" +
			"samConsensus individual.sam ref.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var substitutionsOnly *bool = flag.Bool("substitutionsOnly", false, "This option ignores insertions and deletions and only edits substitutions in the output file.")
	var vcfOutFile *string = flag.String("vcfOutFile", "", "Also write a vcf file of called variants from the reads.")
	var chainOutFile *string = flag.String("chainOutFile", "", "TODO: Also write a chain file describing the positional relationship between the reference and the output fasta.")
	var insertionThreshold *float64 = flag.Float64("insertionThreshold", 0.1, "Requires the number of observations of an insertion relative to read depth required to call an insertion.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	refFile := flag.Arg(1)
	outFile := flag.Arg(2)

	s := Settings {
		SamFileName: inFile,
		RefFile: refFile,
		OutFile: outFile,
		VcfFile: *vcfOutFile,
		SubstitutionsOnly: *substitutionsOnly,
		ChainFile: *chainOutFile,
		InsertionThreshold: *insertionThreshold,
	}

	samConsensus(s)
}
