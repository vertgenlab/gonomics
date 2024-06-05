// Package gaf provides functions for reading, writing, and manipulating GO Annotation File (GAF) format files.
// This package uses the GAF 2.2 file format as specified at the following URL: http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
package gaf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
	"sync"
)

type Gaf struct {
	Db                  string // Database from which DbObjectId is derived. Ex. UniProtKB
	DbObjectId          string // Unique identifier from Db that is associated with DbObjectSymbol
	DbObjectSymbol      string // Gene or protein name.
	Qualifier           string // Relationship between DbObjectSymbol and GoId. Ex. involved_in or located_in. Note that "NOTinvolved_in" implies negation.
	GoId                string // The GO identifier attributed to DbObjectId.
	DbReference         string // Source cited as authority for the attribution of GoId to DbObjectId.
	EvidenceCode        string // Experimental evidence code
	WithFrom            string // Optional: An additional annotation field.
	Aspect              string // One of the following: P (biological process), F (molecular function), or C (cellular component).
	DbObjectName        string // Optional: Name of the gene or gene product. (Usually spelled out, like  T cell receptor gamma variable 8 for TRGV8.
	DbObjectSynonym     string // Optional: additional possible names for the gene, pipe separated. ex. YFL039C|ABY1|END7|actin.
	DbObjectType        string // Type of gene/protein. One of the following: protein_complex; protein; transcript; ncRNA; rRNA; tRNA; snRNA; snoRNA
	Taxon               string // Encodes the species ID.
	Date                string // Date on which the annotation was made. YYYYMMDD format.
	AssignedBy          string // Group who generated annotation. Often matching Db. Differs from Db if made in one database and incorporated in another.
	AnnotationExtension string // Optional. Ex. part_of(CL:0000576)
	GeneProductFormId   string // Optional. Allows annotation of specific splice variants.
}

// Header encodes the header of a Gaf file.
type Header struct {
	Text []string
}

// String converts a Gaf struct into a string for printing.
func (g Gaf) String() string {
	return ToString(g)
}

// ToString converts a Gaf struct into a GAF file format string. Useful for writing to files or printing.
func ToString(g Gaf) string {
	return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
		g.Db,
		g.DbObjectId,
		g.DbObjectSymbol,
		g.Qualifier,
		g.GoId,
		g.DbReference,
		g.EvidenceCode,
		g.WithFrom,
		g.Aspect,
		g.DbObjectName,
		g.DbObjectSynonym,
		g.DbObjectType,
		g.Taxon,
		g.Date,
		g.AssignedBy,
		g.AnnotationExtension,
		g.GeneProductFormId,
	)
}

// ReadHeader processes the contiguous header from an EasyReader
// and advances the Reader past the header lines.
func ReadHeader(file *fileio.EasyReader) Header {
	var answer Header
	var done bool
	var line string
	var err error
	var peek []byte

	for peek, err = file.Peek(1); err == nil && peek[0] == '!' && !done; peek, err = file.Peek(1) {
		line, done = fileio.EasyNextLine(file)
		answer.Text = append(answer.Text, line)
	}

	if err != nil && err != io.EOF {
		log.Fatalf("Error: had the following problem while reading the sam header: %s\n", err)
	}

	return answer
}

// WriteGaf writes an input Gaf struct to an io.Writer.
func WriteGaf(file io.Writer, input Gaf) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", input)
	exception.PanicOnErr(err)
}

// Write writes a slice of Gaf structs to a specified filename.
func Write(filename string, records []Gaf, header Header) {
	var err error
	file := fileio.EasyCreate(filename)

	WriteHeaderToFileHandle(file, header)
	for i := range records {
		WriteGaf(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

// WriteHeaderToFileHandle writes a Gaf header to the input file.
func WriteHeaderToFileHandle(file io.Writer, header Header) {
	var err error
	for i := range header.Text {
		_, err = fmt.Fprintln(file, header.Text[i])
		exception.PanicOnErr(err)
	}
}

// Read returns a slice of Gaf structs from an input filename.
func Read(filename string) ([]Gaf, Header) {
	var line string
	var answer []Gaf
	var err error
	var doneReading bool

	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processGafLine(line)
		answer = append(answer, current)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer, header
}

// processGafLine is a helper function of Read that returns a Gaf struct from an input line of a file.
func processGafLine(line string) Gaf {
	words := strings.Split(line, "\t")
	if len(words) < 15 {
		log.Fatalf("Error: expected at least 15 fields in Gaf line. Found %v on this line:\n%s\n", len(words), line)
	}
	if len(words) > 17 {
		log.Fatalf("Error: expected at most 17 fields in Gaf line. FOund %v opn this line:\n%s\n", len(words), line)
	}
	var answer Gaf = Gaf{
		Db:              words[0],
		DbObjectId:      words[1],
		DbObjectSymbol:  strings.ToUpper(words[2]),
		Qualifier:       words[3],
		GoId:            words[4],
		DbReference:     words[5],
		EvidenceCode:    words[6],
		WithFrom:        words[7],
		Aspect:          words[8],
		DbObjectName:    words[9],
		DbObjectSynonym: words[10],
		DbObjectType:    words[11],
		Taxon:           words[12],
		Date:            words[13],
		AssignedBy:      words[14],
	}
	if len(words) == 16 {
		answer.AnnotationExtension = words[15]
	}
	if len(words) == 17 {
		answer.GeneProductFormId = words[16]
	}
	return answer
}

// NextGaf returns a Gaf struct from an input fileio.EasyReader. REturns a bool that is true when the reader is done.
func NextGaf(reader *fileio.EasyReader) (Gaf, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Gaf{}, true
	}
	return processGafLine(line), false
}

// ReadToChan reads from a fileio.EasyReader to send Gaf structs to a chan<- Gaf.
func ReadToChan(file *fileio.EasyReader, data chan<- Gaf, wg *sync.WaitGroup, header chan<- Header) {
	header <- ReadHeader(file)
	close(header)
	for curr, done := NextGaf(file); !done; curr, done = NextGaf(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GoReadToChan reads Gaf entries from an input filename to a <-chan Gaf.
func GoReadToChan(filename string) (<-chan Gaf, Header) {
	data := make(chan Gaf, 1000)
	header := make(chan Header)

	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	wg.Add(1)
	go ReadToChan(file, data, &wg, header)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, <-header
}
