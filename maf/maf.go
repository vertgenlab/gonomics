// Package maf provides structs and functions for processing alignments in maf format.
// The file format is described more here: https://genome.ucsc.edu/FAQ/FAQformat.html
package maf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"io"
	"log"
	"strings"
	"unicode/utf8"
)

// MafSLine holds the information stored in an "S" line of the maf format.
type MafSLine struct {
	Src     string
	Start   int
	Size    int
	Strand  bool // true for positive strand
	SrcSize int
	Seq     []dna.Base
}

// MafILine holds the information stored in an "I" line of the maf format.
type MafILine struct {
	Src         string
	LeftStatus  rune
	LeftCount   int
	RightStatus rune
	RightCount  int
}

// MafELine holds the information stored in an "E" line of the maf format.
type MafELine struct {
	Src     string
	Start   int
	Size    int
	Strand  bool // true for positive strand
	SrcSize int
	Status  rune
}

// MafSpecies holds the info for a given species in a given maf block
type MafSpecies struct {
	Src   string
	SLine *MafSLine
	ILine *MafILine
	ELine *MafELine
}

// Maf holds the information for a given maf block.
type Maf struct {
	Score   float64
	Species []*MafSpecies
}

// SrcToAssemblyAndChrom many times in the maf format a sequence is labeled
// as assembly.chromosome (e.g. hg38.chr7) and this function breaks
// them apart and returns the assembly followed by the chromosome.
// If there is no dot, it returns the input string as the assembly
// and no chromosome.  If there is more than one dot it calls log.Fatal().
func SrcToAssemblyAndChrom(src string) (string, string) {
	dots := strings.Count(src, ".")
	switch dots {
	case 0:
		return src, ""
	case 1:
		parts := strings.Split(src, ".")
		return parts[0], parts[1]
	default:
		log.Fatalf("Error: too many dots within maf src: %s\n", src)
		return "", ""
	}
}

func parseMafALine(line string) *Maf {
	if !strings.HasPrefix(line, "a") {
		log.Fatalf("Error: the following line should have an 'a' as the first letter: %s\n", line)
	}
	curr := Maf{Score: 0, Species: nil}
	words := strings.Fields(line)
	for i := 1; i < len(words); i++ {
		parts := strings.Split(words[i], "=")
		if parts[0] == "score" {
			curr.Score = parse.StringToFloat64(parts[1])
		}
	}
	return (&curr)
}

func parseMafSLine(line string) *MafSLine {
	words := strings.Fields(line)
	if !strings.HasPrefix(line, "s") || len(words) != 7 {
		log.Fatalf("Error: the following line should start with an 's' and have 7 words: %s\n", line)
	}
	curr := MafSLine{
		Src:     words[1],
		Start:   parse.StringToInt(words[2]),
		Size:    parse.StringToInt(words[3]),
		Strand:  parse.StringToStrand(words[4]),
		SrcSize: parse.StringToInt(words[5]),
		Seq:     dna.StringToBases(words[6]),
	}
	return (&curr)
}

func parseMafIStatus(s string) rune {
	switch s {
	case "C":
		return 'C'
	case "I":
		return 'I'
	case "N":
		return 'N'
	case "n":
		return 'n'
	case "M":
		return 'M'
	case "T":
		return 'T'
	default:
		log.Fatalf("Error: unexpected status for 'i' line in a Maf: %s\n", s)
		return 'X'
	}
}

func parseMafILine(line string) *MafILine {
	words := strings.Fields(line)
	if !strings.HasPrefix(line, "i") || len(words) != 6 {
		log.Fatalf("Error: the following line should start with an 'i' and have 6 words: %s\n", line)
	}
	curr := MafILine{
		Src:         words[1],
		LeftStatus:  parseMafIStatus(words[2]),
		LeftCount:   parse.StringToInt(words[3]),
		RightStatus: parseMafIStatus(words[4]),
		RightCount:  parse.StringToInt(words[5]),
	}
	return (&curr)
}

func parseMafEStatus(s string) rune {
	switch s {
	case "C":
		return 'C'
	case "I":
		return 'I'
	case "M":
		return 'M'
	case "n":
		return 'n'
	case "T": //raven added the "T" case because maf lines can now start with e and end with T
		return 'T'
	default:
		log.Fatalf("Error: unexpected status for 'e' line in a Maf: %s\n", s)
		return 'X'
	}
}

func parseMafELine(line string) *MafELine {
	words := strings.Fields(line)
	if !strings.HasPrefix(line, "e") || len(words) != 7 {
		log.Fatalf("Error: the following line should start with an 'e' and have 6 words: %s\n", line)
	}
	curr := MafELine{
		Src:     words[1],
		Start:   parse.StringToInt(words[2]),
		Size:    parse.StringToInt(words[3]),
		Strand:  parse.StringToStrand(words[4]),
		SrcSize: parse.StringToInt(words[5]),
		Status:  parseMafEStatus(words[6]),
	}
	return (&curr)
}

// FindSpeciesExactMatch take a pointer to a maf block and the string
// to match.  If the source of the maf block is exactly equal to the
// input string than a pointer to the info for that species is returned,
// and otherwise nil is returned.  For example, hg38 would find a maf
// block with hg38 as the source of a sequence, but not hg38.chr7
func FindSpeciesExactMatch(m *Maf, src string) *MafSpecies {
	for i := 0; i < len(m.Species); i++ {
		if m.Species[i].Src == src {
			return m.Species[i]
		}
	}
	return nil
}

// FindSpeciesBeforeDot is similar to FindSpeciesExactMatch, but
// in this case searching hg38 will find blocks with hg38, hg38.chr6
// or hg38.chr22
func FindSpeciesBeforeDot(m *Maf, assembly string) *MafSpecies {
	for i := 0; i < len(m.Species); i++ {
		currAssembly, _ := SrcToAssemblyAndChrom(m.Species[i].Src)
		if currAssembly == assembly {
			return m.Species[i]
		}
	}
	return nil
}

// Read takes the filename of a maf file and returns the contents of the file as a slice of pointers to maf blocks.
// Would be better if this checked for the EOF line at the end and that there was a blank line as the next to last line.
func Read(filename string) []*Maf {
	var answer []*Maf
	var line, prevLine string
	var doneReading bool = false
	var words []string
	var curr *Maf
	var currSpecies *MafSpecies

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, "a") {
			if curr != nil {
				log.Fatalf("Error: no blank line before another 'a' line at line: %s\n", line)
			}
			curr = parseMafALine(line)
		} else if strings.HasPrefix(line, "s") || strings.HasPrefix(line, "i") || strings.HasPrefix(line, "e") {
			if curr == nil {
				log.Fatalf("Error: did not find an 'a' line before this, 'sie' line: %s\n", line)
			}
			words = strings.Fields(line)
			currSpecies = FindSpeciesExactMatch(curr, words[1])
			if currSpecies == nil {
				currSpecies = &MafSpecies{Src: words[1]}
				curr.Species = append(curr.Species, currSpecies)
			}
			if strings.HasPrefix(line, "s") {
				if currSpecies.SLine != nil {
					log.Fatalf("Error: this 's' line looks like a duplicate: %s\n", line)
				}
				currSpecies.SLine = parseMafSLine(line)
			} else if strings.HasPrefix(line, "i") {
				if currSpecies.ILine != nil {
					log.Fatalf("Error: this 'i' line looks like a duplicate: %s\n", line)
				}
				currSpecies.ILine = parseMafILine(line)
			} else if strings.HasPrefix(line, "e") {
				if currSpecies.ELine != nil {
					log.Fatalf("Error: this 'e' line looks like a duplicate: %s\n", line)
				}
				currSpecies.ELine = parseMafELine(line)
			} else {
				log.Fatalf("Error: trouble parsing maf line: %s\n", line)
			}
		} else if line == "" { //blank line at end of maf block
			answer = append(answer, curr)
			curr = nil
		} else {
			log.Fatalf("Unexpected format in maf file on line: %s\n", line)
		}
		prevLine = line
	}
	if prevLine != "" {
		log.Fatalf("Error: maf should have a blank line as the last non-comment line, but found this at end: %s\n", prevLine)
	}
	return answer
}

func calculateFieldSizes(m *Maf) (int, int, int, int) {
	var srcLen, startLen, sizeLen, srcSizeLen int = 1, 1, 1, 1
	var temp int
	for i := 0; i < len(m.Species); i++ {
		if m.Species[i].SLine != nil {
			temp = utf8.RuneCountInString(m.Species[i].SLine.Src)
			if temp > srcLen {
				srcLen = temp
			}
			temp = numbers.DigitsBaseTen(m.Species[i].SLine.Start)
			if temp > startLen {
				startLen = temp
			}
			temp = numbers.DigitsBaseTen(m.Species[i].SLine.Size)
			if temp > sizeLen {
				sizeLen = temp
			}
			temp = numbers.DigitsBaseTen(m.Species[i].SLine.SrcSize)
			if temp > srcSizeLen {
				srcSizeLen = temp
			}
		}
		if m.Species[i].ILine != nil {
			temp = utf8.RuneCountInString(m.Species[i].ILine.Src)
			if temp > srcLen {
				srcLen = temp
			}
		}
		if m.Species[i].ELine != nil {
			temp = utf8.RuneCountInString(m.Species[i].ELine.Src)
			if temp > srcLen {
				srcLen = temp
			}
			temp = numbers.DigitsBaseTen(m.Species[i].ELine.Start)
			if temp > startLen {
				startLen = temp
			}
			temp = numbers.DigitsBaseTen(m.Species[i].ELine.Size)
			if temp > sizeLen {
				sizeLen = temp
			}
			temp = numbers.DigitsBaseTen(m.Species[i].ELine.SrcSize)
			if temp > srcSizeLen {
				srcSizeLen = temp
			}
		}
	}
	return srcLen, startLen, sizeLen, srcSizeLen
}

// WriteToFileHandle writes a single maf block to an io.Writer
func WriteToFileHandle(file io.Writer, m *Maf) {
	_, err := fmt.Fprintf(file, "a score=%.1f\n", m.Score)
	exception.PanicOnErr(err)
	srcChars, startChars, sizeChars, srcSizeChars := calculateFieldSizes(m)
	for i := 0; i < len(m.Species); i++ {
		if m.Species[i].SLine != nil {
			_, err = fmt.Fprintf(file, "s %-*s %*d %*d %c %*d %s\n",
				srcChars, m.Species[i].SLine.Src,
				startChars, m.Species[i].SLine.Start,
				sizeChars, m.Species[i].SLine.Size,
				parse.StrandToRune(m.Species[i].SLine.Strand),
				srcSizeChars, m.Species[i].SLine.SrcSize,
				dna.BasesToString(m.Species[i].SLine.Seq))
			exception.PanicOnErr(err)
		}
		if m.Species[i].ILine != nil {
			//nothing for now
		}
		if m.Species[i].ELine != nil {
			//nothing for now
		}
	}
	_, err = fmt.Fprintf(file, "\n")
	exception.PanicOnErr(err)
}

// Write writes an entire slice of maf blocks to a file specified by the given filename
func Write(filename string, data []*Maf) {
	var err error
	file := fileio.EasyCreate(filename)
	_, err = fmt.Fprint(file, "##maf version=1\n")
	exception.PanicOnErr(err)
	for i := range data {
		WriteToFileHandle(file, data[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}
