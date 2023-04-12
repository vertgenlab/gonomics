// Package fasta provides functions for reading, writing, and manipulating fasta files.
package fasta

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

// Fasta stores the name and sequence of each '>' delimited record in a fasta file.
type Fasta struct {
	Name string
	Seq  []dna.Base
}

// FastaMap stores fasta sequences as a map keyed by the sequence name instead of a slice.
// This allows for easy fasta lookups of chromosomes provided by other files (e.g. BED files).
// A FastaMap can be generated using the ToMap function (e.g. fasta.ToMap(fasta.Read('filename'))).
type FastaMap map[string][]dna.Base

// Read in a fasta file to a []Fasta struct. All sequence records must be preceded by a name line starting with '>'.
// Each record must have a unique sequence name
func Read(filename string) []Fasta {
	return read(filename, NextFasta)
}

// ReadForced functions identically to Read, but any invalid characters in the sequence will be masked to N.
func ReadForced(filename string) []Fasta {
	return read(filename, NextFastaForced)
}

// read is a helper function for Read and ReadForced to reduce duplication
func read(filename string, nextFunc func(*fileio.EasyReader) (Fasta, bool)) []Fasta {
	var curr Fasta
	var answer []Fasta
	var doneReading bool
	usedSeqNames := make(map[string]bool)

	file := fileio.EasyOpen(filename)

	for curr, doneReading = nextFunc(file); !doneReading; curr, doneReading = nextFunc(file) {
		if usedSeqNames[curr.Name] {
			log.Fatalf("ERROR: %s is used as the name for multiple records. Names must be unique.", curr.Name)
		} else {
			usedSeqNames[curr.Name] = true
			answer = append(answer, curr)
		}
	}

	exception.PanicOnErr(file.Close())
	return answer
}

// ReadToString reads a fasta file to a map of sequence strings keyed by the record name.
func ReadToString(filename string) map[string]string {
	ans := make(map[string]string)
	buf := new(bytes.Buffer)
	file := fileio.EasyOpen(filename)

	var line, currName string
	var done bool
	for line, done = fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			if buf.Len() > 0 {
				ans[currName] = buf.String()
				buf.Reset()
			}
			currName = line[1:]
			continue
		}

		buf.WriteString(line)
	}

	if buf.Len() > 0 {
		ans[currName] = buf.String()
	}

	err := file.Close()
	exception.PanicOnErr(err)
	return ans
}

// NextFasta reads a single fasta record from an input EasyReader. Returns true when the file is fully read.
func NextFasta(file *fileio.EasyReader) (Fasta, bool) {
	return nextFasta(file, dna.StringToBases)
}

// NextFastaForced functions identically to Read, but any invalid characters in the sequence will be masked to N.
func NextFastaForced(file *fileio.EasyReader) (Fasta, bool) {
	return nextFasta(file, dna.StringToBasesForced)
}

// nextFasta is a helper function for NextFasta and NextFastaForced to reduce code duplication
func nextFasta(file *fileio.EasyReader, convFunc func(string) []dna.Base) (Fasta, bool) {
	var line string
	var name string
	var answer Fasta
	var doneReading bool
	var peek []byte
	var err error

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			name = line[1:]
			answer = Fasta{Name: name, Seq: make([]dna.Base, 0)}
		} else {
			if answer.Name == "" {
				log.Fatalf("ERROR: %s is missing a sequence name (e.g. >chr1)", file.File.Name())
			}
			answer.Seq = append(answer.Seq, convFunc(line)...)
		}

		peek, err = fileio.EasyPeekReal(file, 1)
		if errors.Is(err, io.EOF) {
			return answer, doneReading // we could return true here, but we do not in order to match the behavior of other Next functions.
		} else {
			exception.PanicOnErr(err)
		}

		if peek[0] == '>' {
			return answer, doneReading
		}
	}
	return answer, doneReading
}

// ToMap converts the a slice of fasta records (e.g. the output of the Read function)
// to a map of sequences keyed to the sequences name.
func ToMap(ref []Fasta) map[string][]dna.Base {
	m := make(map[string][]dna.Base)
	for i := range ref {
		_, ok := m[ref[i].Name]
		if !ok {
			m[ref[i].Name] = ref[i].Seq
		} else {
			log.Panicf("%s used for multiple fasta records. record names must be unique.", ref[i].Name)
		}
	}
	return m
}

// Write a fasta to input filename. Output fastas have line length of 50.
func Write(filename string, records []Fasta) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	WriteToFileHandle(file, records, lineLength)
	exception.PanicOnErr(file.Close())
}

// WriteToFileHandle writes a slice of fasta records to a given io.Writer instead of creating
// a new io.Writer as is done in the Write function.
func WriteToFileHandle(file io.Writer, records []Fasta, lineLength int) {
	for _, rec := range records {
		WriteFasta(file, rec, lineLength)
	}
}

// WriteFasta writes a single fasta record to an io.Writer.
func WriteFasta(file io.Writer, rec Fasta, lineLength int) {
	var err error
	_, err = fmt.Fprintf(file, ">%s\n", rec.Name)
	exception.PanicOnErr(err)
	for i := 0; i < len(rec.Seq); i += lineLength {
		if i+lineLength > len(rec.Seq) {
			_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:]))
			exception.PanicOnErr(err)
		} else {
			_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:i+lineLength]))
			exception.PanicOnErr(err)
		}
	}
}

// Extract will subset a sequence in a fasta file and return a new fasta record with
// the same name and a subset of the sequence. Input start and end are left-closed right-open.
func Extract(f Fasta, start int, end int, name string) Fasta {
	var ans Fasta
	ans.Seq = dna.Extract(f.Seq, start, end)
	ans.Name = name
	return ans
}

// ExtractMulti extracts a subsequence from a fasta file for every entry in a multiFa alignment.
func ExtractMulti(records []Fasta, start int, end int) []Fasta {
	var ans = make([]Fasta, len(records))
	for i := range records {
		ans[i] = Extract(records[i], RefPosToAlnPos(records[0], start), RefPosToAlnPos(records[0], end), records[i].Name)
	}
	return ans
}

// CreateAllGaps creates a fasta record where the sequence is all gaps of length numGaps.
func CreateAllGaps(name string, numGaps int) Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllGaps(numGaps)}
	return answer
}

// CreateAllNs creates a fasta record where the sequence is all Ns of length numN
func CreateAllNs(name string, numN int) Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllNs(numN)}
	return answer
}
