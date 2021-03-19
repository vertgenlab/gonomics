// Package gtf contains functions for reading, writing, and manipulating GTF format files.
// More information on the GTF file format can be found at http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
// Structs in the GTF package are organized hierarchically, with the gene struct containing the underlying transcripts, exons, and other gene features associated with that gene.
package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

// The Gene struct organizes all underlying data on a gene feature in a GTF file.
type Gene struct {
	GeneID      string
	GeneName    string
	Transcripts []*Transcript
}

// The Transcript struct contains information on the location, score and strand of a transcript, along with the underlying exons.
type Transcript struct {
	Chr          string
	Source       string
	Start        int
	End          int
	Score        float64
	Strand       bool
	TranscriptID string
	Exons        []*Exon
}

// The Exon struct contains information on the location, score, and relative order of exons in a GTF file.
type Exon struct {
	Start      int
	End        int
	Score      float64
	ExonNumber string
	ExonID     string
	Cds        *Cds
	FiveUtr    *FiveUtr
	ThreeUtr   *ThreeUtr
}

// FiveUtr contains the location and score information for FiveUtr lines of a GTF file.
type FiveUtr struct {
	Start int
	End   int
	Score float64
}

// Cds contains the location and score information for Cds lines of a GTF file. Cds structs also point to the next and previous Cds in the transcript.
type Cds struct {
	Start int
	End   int
	Score float64
	Frame int
	Prev  *Cds
	Next  *Cds
}

// ThreeUtr contains the location and score information for ThreeUtr lines of a GTF file.
type ThreeUtr struct {
	Start int
	End   int
	Score float64
}

// parseFrame is a helper function for Read that converts a string into the frame value for a Cds struct.
func parseFrame(s string) int {
	if s == "." {
		return -1
	}
	answer := common.StringToInt(s)
	if answer > 2 || answer < 0 {
		log.Fatalf("Frame for GTF entries must be either dot, 0, 1, or 2.")
	}
	return answer
}

// getIds parses identifying lines from a gtfLine
func getIds(words []string) (currGeneID, currGeneName, currT, currEID, currENumber string) {
	att := strings.Split(words[8], ";")
	for i := 0; i < len(att); i++ {
		//this command trims the leading space in the annotation field
		att[i] = strings.TrimSpace(att[i])
		field := strings.Split(att[i], " ")
		if field[0] == "gene_id" {
			currGeneID = strings.Trim(field[1], "\"")
		}
		if field[0] == "transcript_id" {
			currT = strings.Trim(field[1], "\"")
		}
		if field[0] == "gene_name" {
			currGeneName = strings.Trim(field[1], "\"")
		}
		if field[0] == "exon_id" {
			currEID = strings.Trim(field[1], "\"")
		}
		if field[0] == "exon_number" {
			currENumber = strings.Trim(field[1], "\"")
		}
	}
	return
}

// findTranscript finds the transcript with a given ID in a slice of transcripts
func findTranscript(query string, transcripts []*Transcript) *Transcript {
	for i := range transcripts {
		if transcripts[i].TranscriptID == query {
			return transcripts[i]
		}
	}
	return nil
}

// findExon finds the exon witha given ID in a Transcript
func findExon(query string, transcript *Transcript) *Exon {
	for i := range transcript.Exons {
		if transcript.Exons[i].ExonID == query {
			return transcript.Exons[i]
		}
	}
	return nil
}

// parseGtfLine processes a single line of a GTF file and fills in the appropriate data in a Gene struct
func parseGtfLine(line string, currentTranscript *Transcript, prevCds *Cds, answer map[string]*Gene) (*Transcript, *Cds) {
	words := strings.Split(line, "\t")

	if len(words) > 10 {
		log.Fatalf("The GTF file format is limited to nine fields. Line had %d fields.", len(words))
	}

	if words[5] == "." {
		words[5] = "-1"
	}

	currGeneID, currGeneName, currT, currEID, currENumber := getIds(words)

	switch words[2] {
	case "transcript":
		prevCds = nil
		currentTranscript = &Transcript{Chr: words[0], Source: words[1], Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5]), TranscriptID: currT}
		currentTranscript.Strand = common.StringToStrand(words[6])
		currentTranscript.Exons = make([]*Exon, 0)

		if _, ok := answer[currGeneID]; ok {
			answer[currGeneID].Transcripts = append(answer[currGeneID].Transcripts, currentTranscript)
		} else {
			answer[currGeneID] = &Gene{GeneID: currGeneID, GeneName: currGeneName}
			answer[currGeneID].Transcripts = make([]*Transcript, 0)
			answer[currGeneID].Transcripts = append(answer[currGeneID].Transcripts, currentTranscript)
		}
	case "exon":
		t := findTranscript(currT, answer[currGeneID].Transcripts)
		t.Exons = append(t.Exons, &Exon{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), ExonNumber: currENumber, ExonID: currEID, Score: common.StringToFloat64(words[5])})

	case "CDS":
		e := findExon(currEID, findTranscript(currT, answer[currGeneID].Transcripts))
		currentCDS := Cds{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5]), Frame: parseFrame(words[7])}
		currentCDS.Prev = prevCds
		if prevCds != nil {
			prevCds.Next = &currentCDS
		}
		prevCds = &currentCDS
		e.Cds = &currentCDS


	case "5UTR":
		e := findExon(currEID, findTranscript(currT, answer[currGeneID].Transcripts))
		current5Utr := FiveUtr{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5])}
		e.FiveUtr = &current5Utr

	case "3UTR":
		e := findExon(currEID, findTranscript(currT, answer[currGeneID].Transcripts))
		current3Utr := ThreeUtr{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5])}
		e.ThreeUtr = &current3Utr

	default:
		// start_codon and stop_codon lines not read for now.
		// TODO: add in a parser for these lines and throw a log.Fatalf for other line types.
		return currentTranscript, prevCds
	}
	return currentTranscript, prevCds
}

// TODO: Set up Exon and Cds pointers to match the style of transcripts
// Read generates a map[geneID]*Gene of GTF information from an input GTF format file.
func Read(filename string) map[string]*Gene {
	file := fileio.EasyOpen(filename)
	var line string
	var currentTranscript *Transcript
	var doneReading bool
	answer := make(map[string]*Gene)
	var prevCds *Cds

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		currentTranscript, prevCds = parseGtfLine(line, currentTranscript, prevCds, answer)
	}

	err := file.Close()
	exception.PanicOnErr(err)
	return answer
}

// Write writes information contained in a GTF data structure to an output file.
func Write(filename string, records map[string]*Gene) {
	file := fileio.EasyCreate(filename)
	for _, k := range records { //for each gene
		WriteToFileHandle(file, k)
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle is a helper function of Write that facilitates writing GTF data to an output file.
func WriteToFileHandle(file io.Writer, gene *Gene) {
	var err error
	for i := 0; i < len(gene.Transcripts); i++ { //for each transcript associated with that gene
		_, err = fmt.Fprintf(file, "%s\n", gtfTranscriptToString(gene.Transcripts[i], gene))
		exception.PanicOnErr(err)
		for j := 0; j < len(gene.Transcripts[i].Exons); j++ {
			_, err = fmt.Fprintf(file, "%s\n", gtfExonToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
			exception.PanicOnErr(err)
			if gene.Transcripts[i].Exons[j].FiveUtr != nil { //if cds, 5utr, and 3utr are not nil pointers the underlying struct
				_, err = fmt.Fprintf(file, "%s\n", gtf5UtrToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
				exception.PanicOnErr(err)
			}
			if gene.Transcripts[i].Exons[j].Cds != nil {
				_, err = fmt.Fprintf(file, "%s\n", gtfCdsToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
				exception.PanicOnErr(err)
			}
			if gene.Transcripts[i].Exons[j].ThreeUtr != nil {
				_, err = fmt.Fprintf(file, "%s\n", gtf3UtrToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
				exception.PanicOnErr(err)
			}
		}
	}
}

// gtfTranscriptToString is a helper function of WriteToFileHandle that converts a transcript struct into a string to be written to the output file.
func gtfTranscriptToString(t *Transcript, g *Gene) string {
	lineType := "transcript"
	var score, strand, frame, att string
	if t.Score == -1 {
		score = "."
	} else {
		score = fmt.Sprintf("%f", t.Score)
	}
	if t.Strand {
		strand = "+"
	} else {
		strand = "-"
	}
	frame = "."
	att = fmt.Sprintf("gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";", g.GeneID, t.TranscriptID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%s\t%s", t.Chr, t.Source, lineType, t.Start, t.End, score, strand, frame, att)
}

// gtfExonToString is a helper function of WriteToFileHandle that converts an Exon struct into a string to be written to the output file.
func gtfExonToString(e *Exon, t *Transcript, g *Gene) string {
	lineType := "exon"
	var score, strand, frame, att string
	if e.Score == -1 {
		score = "."
	} else {
		score = fmt.Sprintf("%f", e.Score)
	}
	if t.Strand {
		strand = "+"
	} else {
		strand = "-"
	}
	frame = "."
	att = fmt.Sprintf("gene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\"; exon_id \"%s\"; gene_name \"%s\";", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%s\t%s", t.Chr, t.Source, lineType, e.Start, e.End, score, strand, frame, att)
}

// gtf5UtrToString is a helper function of WriteToFileHandle that converts a FiveUtr struct into a string to be written to the output file.
func gtf5UtrToString(e *Exon, t *Transcript, g *Gene) string {
	lineType := "5UTR"
	var score, strand, frame, att string
	if e.FiveUtr.Score == -1 {
		score = "."
	} else {
		score = fmt.Sprintf("%f", e.FiveUtr.Score)
	}
	if t.Strand {
		strand = "+"
	} else {
		strand = "-"
	}
	frame = "."
	att = fmt.Sprintf("gene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\"; exon_id \"%s\"; gene_name \"%s\";", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%s\t%s", t.Chr, t.Source, lineType, e.FiveUtr.Start, e.FiveUtr.End, score, strand, frame, att)
}

// gtfCdsToString is a helper function of WriteToFileHandle that converts a Cds struct into a string to be written to the output file.
func gtfCdsToString(e *Exon, t *Transcript, g *Gene) string {
	lineType := "CDS"
	var score, strand, att string
	if e.Cds.Score == -1 {
		score = "."
	} else {
		score = fmt.Sprintf("%f", e.Cds.Score)
	}
	if t.Strand {
		strand = "+"
	} else {
		strand = "-"
	}
	att = fmt.Sprintf("gene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\"; exon_id \"%s\"; gene_name \"%s\";", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%v\t%s", t.Chr, t.Source, lineType, e.Cds.Start, e.Cds.End, score, strand, e.Cds.Frame, att)
}

// gtf3UtrToString is a helper function of WriteToFileHandle that converts a ThreeUtr struct into a string to be written to the output file.
func gtf3UtrToString(e *Exon, t *Transcript, g *Gene) string {
	lineType := "3UTR"
	var score, strand, frame, att string
	if e.ThreeUtr.Score == -1 {
		score = "."
	} else {
		score = fmt.Sprintf("%f", e.ThreeUtr.Score)
	}
	if t.Strand {
		strand = "+"
	} else {
		strand = "-"
	}
	frame = "."
	att = fmt.Sprintf("gene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\"; exon_id \"%s\"; gene_name \"%s\";", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%v\t%s", t.Chr, t.Source, lineType, e.ThreeUtr.Start, e.ThreeUtr.End, score, strand, frame, att)
}
