//Package gtf contains functions for reading, writing, and manipulating GTF format files.
//More information on the GTF file format can be found at http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
//Structs in the GTF package are organized hierarchically, with the gene struct containing the underlying transcripts, exons, and other gene features associated with that gene.
package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

//The Gene struct organizes all underlying data on a gene feature in a GTF file.
type Gene struct {
	GeneID      string
	GeneName    string
	Transcripts []*Transcript
}

//The Transcript struct contains information on the location, score and strand of a transcript, along with the underlying exons.
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

//The Exon struct contains information on the location, score, and relative order of exons in a GTF file.
type Exon struct {
	Start      int
	End        int
	Score      float64
	ExonNumber string
	ExonID     string
	Cds        *CDS
	FiveUtr    *FiveUTR
	ThreeUtr   *ThreeUTR
}

//FiveUTR contains the location and score information for FiveUTR lines of a GTF file.
type FiveUTR struct {
	Start int
	End   int
	Score float64
}

//CDS contains the location and score information for CDS lines of a GTF file. CDS structs also point to the next and previous CDS in the transcript.
type CDS struct {
	Start int
	End   int
	Score float64
	Frame int
	Prev  *CDS
	Next  *CDS
}

//ThreeUTR contains the location and score information for ThreeUTR lines of a GTF file.
type ThreeUTR struct {
	Start int
	End   int
	Score float64
}

//ParseFram is a helper function for Read that converts a string into the frame value for a CDS struct.
func ParseFrame(s string) int {
	if s == "." {
		return -1
	}
	answer := common.StringToInt(s)
	if answer > 2 || answer < 0 {
		log.Fatalf("Frame for GTF entries must be either dot, 0, 1, or 2.")
	}
	return answer
}

//TODO: Break up into helper functions
//TODO: Set up Exon and CDS pointers to match the style of transcripts
//Read generates a map[geneID]*Gene of GTF information from an input GTF format file.
func Read(filename string) map[string]*Gene {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	var line string
	var currentTranscript *Transcript
	var doneReading bool = false
	answer := make(map[string]*Gene)
	var prevCDS *CDS

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")

		if len(words) > 10 {
			log.Fatalf("The GTF file format is limited to nine fields. Line had %d fields.", len(words))
		}

		if words[5] == "." {
			words[5] = "-1"
		}

		att := strings.Split(words[8], ";")
		var currGeneID, currGeneName, currT, currEID, currENumber string
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

		switch words[2] {
		case "transcript":
			prevCDS = nil
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
			currentExon := Exon{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), ExonNumber: currENumber, ExonID: currEID, Score: common.StringToFloat64(words[5])}
			for i := 0; i < len(answer[currGeneID].Transcripts); i++ {
				if answer[currGeneID].Transcripts[i].TranscriptID == currT {
					answer[currGeneID].Transcripts[i].Exons = append(answer[currGeneID].Transcripts[i].Exons, &currentExon)
				}
			}
		case "CDS":
			currentCDS := CDS{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5]), Frame: ParseFrame(words[7])}
			for i := 0; i < len(answer[currGeneID].Transcripts); i++ {
				if answer[currGeneID].Transcripts[i].TranscriptID == currT {
					for j := 0; j < len(answer[currGeneID].Transcripts[i].Exons); j++ {
						if answer[currGeneID].Transcripts[i].Exons[j].ExonID == currEID {
							currentCDS.Prev = prevCDS
							if prevCDS != nil {
								prevCDS.Next = &currentCDS
							}
							prevCDS = &currentCDS
							answer[currGeneID].Transcripts[i].Exons[j].Cds = &currentCDS
						}
					}
				}
			}
		case "5UTR":
			current5Utr := FiveUTR{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5])}
			for i := 0; i < len(answer[currGeneID].Transcripts); i++ {
				if answer[currGeneID].Transcripts[i].TranscriptID == currT {
					for j := 0; j < len(answer[currGeneID].Transcripts[i].Exons); j++ {
						if answer[currGeneID].Transcripts[i].Exons[j].ExonID == currEID {
							answer[currGeneID].Transcripts[i].Exons[j].FiveUtr = &current5Utr
						}
					}
				}
			}
		case "3UTR":
			current3Utr := ThreeUTR{Start: common.StringToInt(words[3]), End: common.StringToInt(words[4]), Score: common.StringToFloat64(words[5])}
			for i := 0; i < len(answer[currGeneID].Transcripts); i++ {
				for j := 0; j < len(answer[currGeneID].Transcripts[i].Exons); j++ {
					if answer[currGeneID].Transcripts[i].Exons[j].ExonID == currEID {
						answer[currGeneID].Transcripts[i].Exons[j].ThreeUtr = &current3Utr
					}
				}
			}
		default:
			//start_codon and stop_codon lines not read for now.
			//TODO: add in a parser for these lines and throw a log.Fatalf for other line types.
			continue
		}
	}
	return answer
}

//Write writes information contained in a GTF data structure to an output file.
func Write(filename string, records map[string]*Gene) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	for _, k := range records { //for each gene
		WriteToFileHandle(file, k)
	}
}

//WriteToFileHandle is a helper function of Write that facilitates writing GTF data to an output file.
func WriteToFileHandle(file io.Writer, gene *Gene) {
	var err error
	for i := 0; i < len(gene.Transcripts); i++ { //for each transcript associated with that gene
		_, err = fmt.Fprintf(file, "%s\n", GtfTranscriptToString(gene.Transcripts[i], gene))
		common.ExitIfError(err)
		for j := 0; j < len(gene.Transcripts[i].Exons); j++ {
			_, err = fmt.Fprintf(file, "%s\n", GtfExonToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
			common.ExitIfError(err)
			if gene.Transcripts[i].Exons[j].FiveUtr != nil { //if cds, 5utr, and 3utr are not nil pointers the underlying struct
				_, err = fmt.Fprintf(file, "%s\n", Gtf5UtrToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
				common.ExitIfError(err)
			}
			if gene.Transcripts[i].Exons[j].Cds != nil {
				_, err = fmt.Fprintf(file, "%s\n", GtfCdsToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
				common.ExitIfError(err)
			}
			if gene.Transcripts[i].Exons[j].ThreeUtr != nil {
				_, err = fmt.Fprintf(file, "%s\n", Gtf3UtrToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
				common.ExitIfError(err)
			}
		}
	}
}

//GtfTranscriptToString is a helper function of WriteToFileHandle that converts a transcript struct into a string to be written to the output file.
func GtfTranscriptToString(t *Transcript, g *Gene) string {
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

//GtfExonToString is a helper function of WriteToFileHandle that converts an Exon struct into a string to be written to the output file.
func GtfExonToString(e *Exon, t *Transcript, g *Gene) string {
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

//Gtf5UtrToString is a helper function of WriteToFileHandle that converts a FiveUTR struct into a string to be written to the output file.
func Gtf5UtrToString(e *Exon, t *Transcript, g *Gene) string {
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

//GtfCdsToString is a helper function of WriteToFileHandle that converts a CDS struct into a string to be written to the output file.
func GtfCdsToString(e *Exon, t *Transcript, g *Gene) string {
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

//Gtf3UtrToString is a helper function of WriteToFileHandle that converts a ThreeUTR struct into a string to be written to the output file.
func Gtf3UtrToString(e *Exon, t *Transcript, g *Gene) string {
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
