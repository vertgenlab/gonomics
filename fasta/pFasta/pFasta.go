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

	WriteHeader(out, records)
	var base pDna.Float32Base
	for _, record := range records {
		for _, base = range record.Seq {
			WriteBase(out, base)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func WriteHeader(out *fileio.EasyWriter, records []PFasta) {
	_, err := fmt.Fprint(out, "pFasta_format_1.0\n")
	exception.PanicOnErr(err)
	for i := range records {
		_, err = fmt.Fprintf(out, "%s\t%v\n", records[i].Name, len(records[i].Seq))
		exception.PanicOnErr(err)
	}
	_, err = fmt.Fprintf(out, "EndHeader\n")
	exception.PanicOnErr(err)
}

func WriteBase(out *fileio.EasyWriter, base pDna.Float32Base) {
	writeBaseProbability(out, base.A)
	writeBaseProbability(out, base.C)
	writeBaseProbability(out, base.G)
	writeBaseProbability(out, base.T)
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

// ToMap converts a slice of pFasta records (e.g. the output of the Read function)
// to a map of sequences keyed to the sequences name.
func ToMap(ref []PFasta) map[string][]pDna.Float32Base {
	m := make(map[string][]pDna.Float32Base)
	for i := range ref {
		_, ok := m[ref[i].Name]
		if !ok {
			m[ref[i].Name] = ref[i].Seq
		} else {
			log.Panicf("%s used for multiple pFasta records. record names must be unique.", ref[i].Name)
		}
	}
	return m
}

/*
// QC tests if each base in a single pFasta sequence is valid. Valid means either gap {0 0 0 0} or non-gap (bases' probabilities add up to 1)
func QC(pfa []pDna.Float32Base) {
	for i := range pfa {
		if !pDna.IsGap(pfa[i]) && !pDna.IsNonGap(pfa[i]) {
			log.Fatalf("Error: pFasta QC failed. Position %v has invalid base: %v\n", i, pfa[i])
		}
	}
}

// MakeValid makes each base in a single pFasta sequence valid
func MakeValid(pfaIn []pDna.Float32Base) []pDna.Float32Base {
	pfaOut := make([]pDna.Float32Base, len(pfaIn))
	for i := range pfaIn {
		if pDna.IsGap(pfaIn[i]) {
			pfaOut[i] = pfaIn[i]
		} else if pDna.IsNonGap(pfaIn[i]) {
			pfaOut[i] = pfaIn[i]
		} else { // already checked that the base is not gap, so convert it to valid non-gap
			pfaOut[i] = pDna.MakeValid(pfaIn[i])
		}
	}
	return pfaOut
}
*/
