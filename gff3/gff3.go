// Package gff3 provides utilities to process GFF3 files.
package gff3

import (
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

// Gff3 represents a single entry in a GFF3 file.
type Gff3 struct {
	Id         string  // Unique identifier for the feature
	Source     string  // Source of the feature (e.g., program or database)
	Type       string  // Type of feature (from SOFA ontology)
	Start      int     // Start position of the feature (1-based)
	End        int     // End position of the feature (1-based)
	Score      float64 // Score or confidence value
	Strand     byte    // Strand (+ or -)
	Phase      string  // Phase for coding sequence features (0, 1, or 2)
	Attributes []Tag   // Additional attributes as tag-value pairs
}

// Tag represents a single tag-value pair in the GFF3 attributes field.
type Tag struct {
	Label string // Tag label (e.g., ID, Name, Parent)
	Value string // Tag value
}

// Read reads a GFF3 file and returns a slice of Gff3 structs.
func Read(filename string) []Gff3 {
	var answer []Gff3
	file := fileio.NewByteReader(filename)

	line, done := fileio.ReadLine(file)

	if !strings.HasPrefix(line.String(), "##gff-version 3") || done {
		log.Panicf("Error: Header line does not contain the text '##gff-version 3' or doesn't contain any data.")
	}

	for line, done = fileio.ReadLine(file); !done; line, done = fileio.ReadLine(file) {
		if line.String() != "###" {
			answer = append(answer, ToGff3(line.String()))
		}
	}
	exception.PanicOnErr(file.Close())
	return answer
}

// Write writes a slice of Gff3 structs with a specified number of fields to a specified filename.
func Write(filename string, records []Gff3) {
	file := fileio.EasyCreate(filename)
	WriteToFileHandle(file, records)
	err := file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes an input Gff3 struct with a specified number of fields to an io.Writer.
func WriteToFileHandle(file io.Writer, records []Gff3) {
	_, err := fmt.Fprintln(file, "##gff-version 3")
	exception.PanicOnErr(err)
	for _, rec := range records {
		_, err = fmt.Fprintf(file, "%s\n", Gff3ToString(rec))
		exception.PanicOnErr(err)
	}
	_, err = fmt.Fprintln(file, "###")
	exception.PanicOnErr(err)
}

// ToGff3 converts a single line from a GFF3 file to a Gff3 struct.
func ToGff3(line string) Gff3 {
	row := strings.Split(line, "\t")
	var err error

	if len(row) != 9 {
		log.Panicf("Error: invalid GFF3 line: expected 9 fields, got %d\n%s", len(row), line)
	}

	ans := Gff3{
		Id:     row[0],
		Source: row[1],
		Type:   row[2],
		Phase:  row[7],
	}

	ans.Start, err = strconv.Atoi(row[3])
	exception.PanicOnErr(err)

	ans.End, err = strconv.Atoi(row[4])
	exception.PanicOnErr(err)

	if row[5] != "." {
		ans.Score, err = strconv.ParseFloat(row[5], 64)
		exception.PanicOnErr(err)
	}

	switch row[6] {
	case "+":
		ans.Strand = '+'
	case "-":
		ans.Strand = '-'
	default:
		log.Panicf("Error: Field 6 does not equal '+' or '-'.")
	}
	ans.Attributes = parseTags(row[8])
	return ans
}

// Gff3ToString converts a Gff3 struct to a string in GFF3 format.
func Gff3ToString(g Gff3) string {
	var buf strings.Builder
	var err error

	_, err = buf.WriteString(g.Id)
	exception.PanicOnErr(err)
	exception.PanicOnErr(buf.WriteByte('\t'))
	_, err = buf.WriteString(g.Source)
	exception.PanicOnErr(err)
	exception.PanicOnErr(buf.WriteByte('\t'))
	_, err = buf.WriteString(g.Type)
	exception.PanicOnErr(err)
	exception.PanicOnErr(buf.WriteByte('\t'))
	_, err = buf.WriteString(strconv.Itoa(g.Start))
	exception.PanicOnErr(err)
	exception.PanicOnErr(buf.WriteByte('\t'))
	_, err = buf.WriteString(strconv.Itoa(g.End))
	exception.PanicOnErr(err)
	exception.PanicOnErr(buf.WriteByte('\t'))

	// Handle potential '.' for score
	if g.Score == 0 {
		exception.PanicOnErr(buf.WriteByte('.'))
	} else {
		_, err = buf.WriteString(fmt.Sprintf("%g", g.Score)) // Use %g for shortest representation
		exception.PanicOnErr(err)
	}
	exception.PanicOnErr(buf.WriteByte('\t'))

	exception.PanicOnErr(buf.WriteByte(g.Strand))
	exception.PanicOnErr(buf.WriteByte('\t'))
	_, err = buf.WriteString(g.Phase)
	exception.PanicOnErr(err)
	exception.PanicOnErr(buf.WriteByte('\t'))

	_, err = buf.WriteString(tagsToString(g.Attributes))
	exception.PanicOnErr(err)

	return buf.String()
}

// parseTags parses a semicolon-separated string of tag-value pairs and returns a slice of Tag structs.
func parseTags(attrStr string) []Tag {
	attributes := []Tag{}
	for _, pairStr := range strings.Split(attrStr, ";") {
		if pairStr == "" {
			continue // Skip empty attributes
		}
		parts := strings.Split(pairStr, "=")
		if len(parts) != 2 {
			log.Panicf("Error: Invalid attribute pair: %s", pairStr)
		}
		attributes = append(attributes, Tag{Label: parts[0], Value: parts[1]})
	}
	return attributes
}

// tagsToString converts attributes from a slice of Tags to string
func tagsToString(attr []Tag) string {
	var buf strings.Builder
	var err error
	for i, tag := range attr {
		_, err = buf.WriteString(tag.Label)
		exception.PanicOnErr(err)
		exception.PanicOnErr(buf.WriteByte('='))
		_, err = buf.WriteString(tag.Value)
		exception.PanicOnErr(err)
		if i < len(attr)-1 {
			exception.PanicOnErr(buf.WriteByte(';'))
		}
	}
	return buf.String()
}
