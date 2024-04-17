package pFasta

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/x448/float16"
	"io"
	"log"
	"strings"
)

// PFasta is the probabilistic analog of fasta, which represents sequences of
// pDna bases rather than dna bases.
type PFasta struct {
	Name string
	Seq  []pDna.Float32Base
}

// Write takes an input []pFasta and writes to an output binary pFa file.
func Write(outFile string, records []PFasta) {
	var err error
	out := fileio.EasyCreate(outFile)

	// first we write the header
	_, err = fmt.Fprint(out, "pFasta_format_1.0\n")
	exception.PanicOnErr(err)
	for i := range records {
		_, err = fmt.Fprintf(out, "%s\t%v\n", records[i].Name, len(records[i].Seq))
		exception.PanicOnErr(err)
	}
	_, err = fmt.Fprintf(out, "EndHeader\n")
	exception.PanicOnErr(err)

	var currSeqIndex, currPos int
	for currSeqIndex = range records {
		for currPos = range records[currSeqIndex].Seq {
			writeBaseProbability(out, records[currSeqIndex].Seq[currPos].A)
			writeBaseProbability(out, records[currSeqIndex].Seq[currPos].C)
			writeBaseProbability(out, records[currSeqIndex].Seq[currPos].G)
			writeBaseProbability(out, records[currSeqIndex].Seq[currPos].T)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

// writeBase is a helper function of Write, which write an input float32
// base probability into an output file as binary.
func writeBaseProbability(out *fileio.EasyWriter, b float32) {
	err := binary.Write(out, binary.LittleEndian, float16.Fromfloat32(b).Bits())
	exception.PanicOnErr(err)
}

// makeEmptyRecords is a helper function of Read. It parses the metadata header of
// a binary pFa file and returns the skeleton of the []PFasta, where the names and sequences
// have been initialized but all probabilities are 0.
func makeEmptyRecords(reader *bufio.Reader) []PFasta {
	var records = make([]PFasta, 0)
	unparsedHeader := ReadPfaHeader(reader)

	if unparsedHeader[0] != "pFasta_format_1.0" {
		log.Fatalf("Error: unrecognized pFasta file format. Found: %v.\n", unparsedHeader[0])
	}

	if unparsedHeader[len(unparsedHeader)-1] != "EndHeader" {
		log.Fatalf("Error: unrecognized end of header. Found: %v.\n", unparsedHeader[len(unparsedHeader)-1])
	}

	var words []string
	for currHeaderLine := 1; currHeaderLine < len(unparsedHeader)-1; currHeaderLine++ {
		words = strings.Split(unparsedHeader[currHeaderLine], "\t")
		records = append(records, PFasta{
			Name: words[0],
			Seq:  make([]pDna.Float32Base, parse.StringToInt(words[1]))})
	}
	return records
}

// ReadPfaHeader parses the header of a pFasta format file from an input *bufio.Reader.
func ReadPfaHeader(reader *bufio.Reader) []string {
	var header []string
	var line string
	var err error

	for { // this loop will end either at the end of the file (with panic), or when 'EndHeader' is reached.
		line, err = reader.ReadString('\n')
		exception.PanicOnErr(err)
		line = strings.TrimSuffix(line, "\n")
		line = strings.TrimSuffix(line, "\r")
		if line == "EndHeader" {
			header = append(header, line)
			break
		}
		header = append(header, line)
	}

	return header
}

// Read parses a []PFasta from an input file handle in .pFa binary format.
func Read(inFile string) []PFasta {
	var err error
	file := fileio.EasyOpen(inFile)
	records := makeEmptyRecords(file.BuffReader)
	var currBase = make([]byte, 2) //8 bytes, or 64 bit slice. Enough for one pDna base probability.
	var currSeq, currPos int

	for currSeq = range records {
		for currPos = range records[currSeq].Seq {
			_, err = io.ReadFull(file.BuffReader, currBase)
			exception.PanicOnErr(err)
			records[currSeq].Seq[currPos].A = float16.Frombits(binary.LittleEndian.Uint16(currBase)).Float32()
			_, err = io.ReadFull(file.BuffReader, currBase)
			exception.PanicOnErr(err)
			records[currSeq].Seq[currPos].C = float16.Frombits(binary.LittleEndian.Uint16(currBase)).Float32()
			_, err = io.ReadFull(file.BuffReader, currBase)
			exception.PanicOnErr(err)
			records[currSeq].Seq[currPos].G = float16.Frombits(binary.LittleEndian.Uint16(currBase)).Float32()
			_, err = io.ReadFull(file.BuffReader, currBase)
			exception.PanicOnErr(err)
			records[currSeq].Seq[currPos].T = float16.Frombits(binary.LittleEndian.Uint16(currBase)).Float32()
		}
	}
	// check to make sure we are at the end of the file
	currBase = currBase[0:1]
	_, err = io.ReadFull(file.BuffReader, currBase)
	if !errors.Is(err, io.EOF) {
		log.Fatalf("Error: %s has more sequence than expected from the header\n", inFile)
	}

	err = file.Close()
	exception.PanicOnErr(err)
	return records
}
