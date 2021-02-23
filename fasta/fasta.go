// Package fasta provides functions for reading, writing, and manipulating fasta files.
package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

// Fasta stores the name and sequence of each '>' delimited record in a fasta file.
type Fasta struct {
	Name string
	Seq  []dna.Base
}

// Read in a fasta file to a []*Fasta struct. All sequence records must be preceded by a name line starting with '>'.
// Each record must have a unique sequence name
func Read(filename string) []*Fasta {
	var line string
	var name string
	var answer []*Fasta
	var seqIdx int = -1
	var doneReading bool = false
	usedSeqNames := make(map[string]bool)

	file := fileio.EasyOpen(filename)

	prefix, err := fileio.EasyPeekReal(file, 1)
	if err != nil || prefix[0] != '>' {
		log.Fatalf("ERROR: %s is missing a sequence name (e.g. >chr1)", filename)
	}

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			name = line[1:]
			if usedSeqNames[name] {
				log.Fatalf("ERROR: %s is used as the name for multiple records. Names must be unique.", name)
			} else {
				usedSeqNames[name] = true
				answer = append(answer, &Fasta{Name: name})
				seqIdx++
			}
		} else {
			answer[seqIdx].Seq = append(answer[seqIdx].Seq, dna.StringToBases(line)...)
		}
	}

	exception.WarningOnErr(file.Close())
	return answer
}

// FastaMap converts the output of the Read function to a map of sequences keyed to the sequences name.
func FastaMap(ref []*Fasta) map[string][]dna.Base {
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
func Write(filename string, records []*Fasta) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	WriteToFileHandle(file, records, lineLength)
	exception.WarningOnErr(file.Close())
}

// WriteToFileHanle writes a slice of fasta records to a given io.Writer instead of creating
// a new io.Writer as is done in the Write function.
func WriteToFileHandle(file io.Writer, records []*Fasta, lineLength int) {
	for _, rec := range records {
		WriteFasta(file, rec, lineLength)
	}
}

// WriteFasta writes a single fasta record to an io.Writer.
func WriteFasta(file io.Writer, rec *Fasta, lineLength int) {
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

// SplitWrite writes each fasta record to its own file.
func SplitWrite(filename string, records []*Fasta) {
	for _, rec := range records {
		Write(filename+"_"+rec.Name+"_"+".fasta", []*Fasta{rec})
	}
}

// CreateAllGaps creates a fasta record where the sequence is all gaps of length numGaps.
func CreateAllGaps(name string, numGaps int) *Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllGaps(numGaps)}
	return &answer
}

// CreateAllNs creates a fasta record where the sequence is all Ns of length numN
func CreateAllNs(name string, numN int) *Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllNs(numN)}
	return &answer
}
