package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

type Gene struct {
	GeneID      string
	GeneName    string
	Transcripts []*Transcript
}

//TODO: Functionality for canonical exon as longest transcript.
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

type FiveUTR struct {
	Start int
	End   int
	Score float64
}

type CDS struct {
	Start int
	End   int
	Score float64
	Frame int
	Prev  *CDS
	Next  *CDS
}

type ThreeUTR struct {
	Start int
	End   int
	Score float64
}

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

//reads to  map[geneID]*Gene
//TODO: Break up into helper functions
//TODO: Set up Exon and CDS pointers to match the style of transcripts
func Read(filename string) map[string]*Gene {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	var line string
	var currentTranscript *Transcript
	var doneReading bool = false
	answer := make(map[string]*Gene)

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
				currGeneID = field[1]
			}
			if field[0] == "transcript_id" {
				currT = field[1]
			}
			if field[0] == "gene_name" {
				currGeneName = field[1]
			}
			if field[0] == "exon_id" {
				currEID = field[1]
			}
			if field[0] == "exon_number" {
				currENumber = field[1]
			}
		}
		var prevCDS *CDS
		switch words[2] {
		case "transcript":
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

func Write(filename string, records map[string]*Gene) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	var err error

	for _, k := range records { //for each gene
		err = WriteToFileHandle(file, k)
		common.ExitIfError(err)
	}
}

func WriteToFileHandle(file io.Writer, gene *Gene) error {
	var err error
	for i := 0; i < len(gene.Transcripts); i++ { //for each transcript associated with that gene
		_, err = fmt.Fprintf(file, "%s\n", GtfTranscriptToString(gene.Transcripts[i], gene))
		for j := 0; j < len(gene.Transcripts[i].Exons); j++ {
			_, err = fmt.Fprintf(file, "%s\n", GtfExonToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
			if gene.Transcripts[i].Exons[j].FiveUtr != nil { //if cds, 5utr, and 3utr are not nil pointers the underlying struct
				_, err = fmt.Fprintf(file, "%s\n", Gtf5UtrToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
			}
			if gene.Transcripts[i].Exons[j].Cds != nil {
				_, err = fmt.Fprintf(file, "%s\n", GtfCdsToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
			}
			if gene.Transcripts[i].Exons[j].ThreeUtr != nil {
				_, err = fmt.Fprintf(file, "%s\n", Gtf3UtrToString(gene.Transcripts[i].Exons[j], gene.Transcripts[i], gene))
			}
		}
	}
	return err
}

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
	att = fmt.Sprintf("gene_id %s; transcript_id %s; gene_name %s;", g.GeneID, t.TranscriptID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%s\t%s", t.Chr, t.Source, lineType, t.Start, t.End, score, strand, frame, att)
}

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
	att = fmt.Sprintf("gene_id %s; transcript_id %s; exon_number %s; exon_id %s; gene_name %s;", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%s\t%s", t.Chr, t.Source, lineType, e.Start, e.End, score, strand, frame, att)
}

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
	att = fmt.Sprintf("gene_id %s; transcript_id %s; exon_number %s; exon_id %s; gene_name %s;", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%s\t%s", t.Chr, t.Source, lineType, e.FiveUtr.Start, e.FiveUtr.End, score, strand, frame, att)
}

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
	att = fmt.Sprintf("gene_id %s; transcript_id %s; exon_number %s; exon_id %s; gene_name %s;", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%v\t%s", t.Chr, t.Source, lineType, e.Cds.Start, e.Cds.End, score, strand, e.Cds.Frame, att)
}

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
	att = fmt.Sprintf("gene_id %s; transcript_id %s; exon_number %s; exon_id %s; gene_name %s;", g.GeneID, t.TranscriptID, e.ExonNumber, e.ExonID, g.GeneName)
	return fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%s\t%s\t%v\t%s", t.Chr, t.Source, lineType, e.ThreeUtr.Start, e.ThreeUtr.End, score, strand, frame, att)
}
