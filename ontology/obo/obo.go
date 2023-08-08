// Package obo provides functions for reading, writing, and manipulating Open Biomedical Ontologies (obo) format files.
// This package follows the OBO Flat File Format 1.4, as specified here: http://owlcollab.github.io/oboformat/doc/obo-syntax.html
package obo

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

// Obo is a struct representing one node in an Obo format file, or one ontology term.
// Required and relevant fields are parsed. As this format is quite flexible and permissive, all
// other fields will be added to 'OtherFields'.
// Go:0000001 from the GO Catalog is used for examples of each field.
type Obo struct {
	Id            string              // Required, can only appear once. GO ID. Ex. GO:0000001
	Name          string              // Required, can only appear once. GO Name. Ex. http://owlcollab.github.io/oboformat/doc/obo-syntax.html
	NameSpace     string              // Required, can only appear once. The category or domain of the term. Ex. biological_process
	Def           string              // Required, can only appear once. Sentence definition of the term.
	IsObsolete    bool                // Optional. Can only appear once. Specifies obsolete ontology labels, which can be filtered downstream.
	IsA           []IsADescription    // Optional. Can appear multiple times in an entry. Specifies parent nodes. Ex. is_a: GO:0048308 ! organelle inheritance // is_a: GO:0048311 ! mitochondrion distribution
	Synonyms      []string            // Optional. Can appear multiple times in an entry. Specifies alternative names for term.
	XRefs         []string            // Optional. Can appear multiple times in an entry. Specifies an external reference about the term. Ex. Wikipedia:Reproduction
	AltIds        []string            // Optional. Can appear multiple times in an entry. Specify an alternative GO ID that also points to this entry.
	Relationships []string            // Optional. Can appear multiple times in an entry. Specifies relationships between ontologies other than 'is_a'.
	Comments      []string            // Optional. Can appear multiple times in an entry. Stores additional comments about an ontology.
	OtherFields   map[string][]string // Catch all for other field names.
	Parents       []*Obo              // Parent nodes of this ontology.
	Children      []*Obo              // Child nodes of this ontology.
	SubTreeSize   int                 //Total number of descendents.
}

// Header encodes the header of an Obo file. Header lines are not marked by a leading character,
// but are defined as lines preceding the first instance of '[Term]'.
type Header struct {
	Text []string
}

// IsADescription contains the information in the "IsA" field of an Obo record.
type IsADescription struct {
	ParentId   string
	ParentInfo []string
}

// String is a method for IsADescription which formats an input struct as a string.
func (d IsADescription) String() string {
	return IsADescriptionToString(d)
}

// IsADescriptionToString converts an IsADescription struct into a string.
func IsADescriptionToString(d IsADescription) string {
	var answer = ""
	answer = answer + d.ParentId
	for i := range d.ParentInfo {
		answer = answer + " " + d.ParentInfo[i]
	}
	return answer
}

// ReadHeader processes the contiguous header from an EasyReader
// and advances the Reader past the header lines. Note that Header lines are not marked by a leading\
// character, but are defined as lines preceding the first instance of '[Term]'.
func ReadHeader(file *fileio.EasyReader) Header {
	var answer Header
	var done bool
	var line string
	var err error
	var peek []byte

	// we will add text to the header as long as we do not encounter the string "[Term]", which signifies a new entry.
	for peek, err = file.Peek(6); err == nil && !peekIsTermLine(peek) && !done; peek, err = file.Peek(6) {
		line, done = fileio.EasyNextLine(file)
		answer.Text = append(answer.Text, line)
	}

	return answer
}

// peekIsTermLine returns true if an input slice of bytes is equal to "[Term]", which is the signifier of a
// new ontology element.
func peekIsTermLine(peek []byte) bool {
	if len(peek) != 6 {
		log.Fatalf("Error: peekIsTermLine expected 6 bytes. Received %v.\n", len(peek))
	}
	peekString := string(peek)
	return peekString == "[Term]"
}

// Read creates a slice of Obo structs and a Header struct from an input filename.
func Read(filename string) ([]Obo, Header) {
	var curr Obo
	var done bool
	var answer []Obo = make([]Obo, 0)

	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	for curr, done = NextObo(file); !done; curr, done = NextObo(file) {
		answer = append(answer, curr)
	}
	return answer, header
}

// NextObo is a helper function of Read and GoReadToChan. NextObo checks a reader for additional lines of data
// and returns an Obo struct and a bool telling the caller if it is done reading from the file.
func NextObo(reader *fileio.EasyReader) (Obo, bool) {
	var currLine string
	var done, endOfEntry bool = false, false
	var lines []string = make([]string, 0)
	for !done && !endOfEntry {
		currLine, done = fileio.EasyNextRealLine(reader)
		if currLine == "" {
			if len(lines) > 0 && lines[0] == "[Typedef]" {
				lines = lines[:0] //clear the slice, we ignore Typedef entries.
			} else {
				endOfEntry = true
			}
		} else {
			lines = append(lines, currLine)
		}
	}
	if done {
		return Obo{}, true
	}
	return processOboTerm(lines), false
}

// processOboTerm is a helper function of NextObo that parses an Obo struct from a set of input lines.
func processOboTerm(lines []string) Obo {
	var words, isAWords []string
	var isAStruct IsADescription
	var answer Obo
	answer.OtherFields = make(map[string][]string)
	for i := range lines {
		words = strings.Split(lines[i], ": ") // sep is colon-space.
		if len(words) == 0 {                  // if we've hit a blank line, we return the current struct.
			return answer
		}
		if len(words) == 1 {
			if words[0] != "[Term]" {
				log.Fatalf("Error: Unrecognized line in the following entry:\n%v", lines)
			}
		} else {
			switch words[0] {
			case "id":
				if answer.Id != "" {
					log.Fatalf("Error: more than one ID found for the following Obo entry:\n%v", lines)
				}
				answer.Id = words[1]
			case "name":
				if answer.Name != "" {
					log.Fatalf("Error: more than one name found for the following Obo entry:\n%v", lines)
				}
				answer.Name = words[1]
			case "namespace":
				if answer.NameSpace != "" {
					log.Fatalf("Error: more than one namespace found for the following Obo entry:\n%v", lines)
				}
				answer.NameSpace = words[1]
			case "def":
				if answer.Def != "" {
					log.Fatalf("Error: more than one def found for the following Obo entry:\n%v", lines)
				}
				answer.Def = words[1]
			case "is_obsolete":
				if words[1] == "true" {
					answer.IsObsolete = true
				} else {
					log.Fatalf("Error: Unrecognized entry after is_obsolete in this Obo entry:\n%v", lines)
				}
			case "is_a":
				isAWords = strings.Split(words[1], " ")
				isAStruct = IsADescription{
					ParentId:   isAWords[0],
					ParentInfo: isAWords[1:],
				}
				answer.IsA = append(answer.IsA, isAStruct)
			case "synonym":
				answer.Synonyms = append(answer.Synonyms, words[1])
			case "xref":
				answer.XRefs = append(answer.XRefs, words[1])
			case "alt_id":
				answer.AltIds = append(answer.AltIds, words[1])
			case "relationship":
				answer.Relationships = append(answer.Relationships, words[1])
			case "comment":
				answer.Comments = append(answer.Comments, words[1])
			default:
				answer.OtherFields[words[0]] = append(answer.OtherFields[words[0]], words[1])
			}
		}
	}

	if answer.Id == "" {
		log.Fatalf("Error: 'id' not found in the following term:\n%v", lines)
	}
	if answer.Name == "" {
		log.Fatalf("Error: 'name' not found in the following term:\n%v", lines)
	}
	if answer.NameSpace == "" {
		log.Fatalf("Error: 'namespace' not found in the following term:\n%v", lines)
	}
	if answer.Def == "" {
		log.Fatalf("Error: 'def' not found in the following term:\n%v", lines)
	}

	return answer
}

func (o Obo) String() string {
	return ToString(o)
}

// ToString converts an Obo struct into an OBO file format string. Used for writing to files or printing.
func ToString(o Obo) string {
	var i int
	var currKey string
	var answer = fmt.Sprintf("[Term]\nid: %s\nname: %s\nnamespace: %s\ndef: %s\n", o.Id, o.Name, o.NameSpace, o.Def)
	if o.IsObsolete {
		answer = answer + "is_obsolete: true\n"
	}
	for i = range o.IsA {
		answer = answer + fmt.Sprintf("is_a: %s\n", o.IsA[i])
	}
	for i = range o.Synonyms {
		answer = answer + fmt.Sprintf("synonym: %s\n", o.Synonyms[i])
	}
	for i = range o.XRefs {
		answer = answer + fmt.Sprintf("xref: %s\n", o.XRefs[i])
	}
	for i = range o.AltIds {
		answer = answer + fmt.Sprintf("alt_id: %s\n", o.AltIds[i])
	}
	for i = range o.Relationships {
		answer = answer + fmt.Sprintf("relationship: %s\n", o.Relationships[i])
	}
	for i = range o.Comments {
		answer = answer + fmt.Sprintf("comment: %s\n", o.Comments[i])
	}
	for currKey = range o.OtherFields {
		for i = range o.OtherFields[currKey] {
			answer = answer + fmt.Sprintf("%s: %s\n", currKey, o.OtherFields[currKey][i])
		}
	}
	return answer
}

// WriteObo writes an input Bed struct to an io.Writer.
func WriteObo(file io.Writer, input Obo) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", input)
	exception.PanicOnErr(err)
}

// WriteHeaderToFileHandle writes an Obo header to an io.Writer.
func WriteHeaderToFileHandle(file io.Writer, header Header) {
	var err error
	for i := range header.Text {
		_, err = fmt.Fprintln(file, header.Text[i])
		exception.PanicOnErr(err)
	}
}

// Write writes a slice of Obo structs to a specified filename.
func Write(filename string, records []Obo, header Header) {
	var err error
	file := fileio.EasyCreate(filename)
	WriteHeaderToFileHandle(file, header)
	for i := range records {
		WriteObo(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}
