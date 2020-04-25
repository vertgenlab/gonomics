package bed

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

type Bed struct {
	Chrom      string
	ChromStart int64
	ChromEnd   int64
	Name       string
	Score      int64
	Strand     bool
	Annotation []string //long form for extra fields
}

func BedToString(bunk *Bed, fields int) string {
	switch fields {
	case 3:
		return fmt.Sprintf("%s\t%d\t%d", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd)
	case 4:
		return fmt.Sprintf("%s\t%d\t%d\t%s", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name)
	case 5:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score)
	case 6:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%c", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score, common.StrandToRune(bunk.Strand))
	case 7:
		var out string = fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%c", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score, common.StrandToRune(bunk.Strand))
		for i := 0; i < len(bunk.Annotation); i++ {
			out = fmt.Sprintf("%s\t%s", out, bunk.Annotation[i])
		}
		return out
	default:
		log.Fatalf("Error: expecting a request to print 3 to 7 bed fields, but got: %d\n", fields)
	}
	return ""
}

func WriteToFileHandle(file io.Writer, rec *Bed, fields int) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", BedToString(rec, fields))
	common.ExitIfError(err)
}

func WriteSliceToFileHandle(file io.Writer, records []*Bed, fields int) {
	for _, rec := range records {
		WriteToFileHandle(file, rec, fields)
	}
}

func Write(filename string, records []*Bed, fields int) {
	file := fileio.EasyCreate(filename)
	defer file.Close()

	WriteSliceToFileHandle(file, records, fields)
}

func Read(filename string) []*Bed {
	var line string
	var answer []*Bed
	var doneReading bool = false

	file := fileio.EasyOpen(filename)
	defer file.Close()
	//reader := bufio.NewReader(file)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processBedLine(line)
		answer = append(answer, current)
	}
	return answer
}

func processBedLine(line string) *Bed {
	words := strings.Split(line, "\t")
	startNum := common.StringToInt64(words[1])
	endNum := common.StringToInt64(words[2])

	current := Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
	if len(words) >= 4 {
		current.Name = words[3]
	}
	if len(words) >= 5 {
		current.Score = common.StringToInt64(words[4])
	}
	if len(words) >= 6 {
		current.Strand = common.StringToStrand(words[5])
	}
	if len(words) >= 7 {
		for i := 6; i < len(words); i++ {
			current.Annotation = append(current.Annotation, words[i])
		}
	}
	return &current
}

func NextBed(reader *fileio.EasyReader) (*Bed, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	return processBedLine(line), false
}

func GoReadToChan(filename string) <-chan *Bed {
	output := make(chan *Bed)
	go ReadToChan(filename, output)
	return output
}

func ReadToChan(filename string, output chan<- *Bed) {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	var curr *Bed
	var done bool
	for curr, done = NextBed(file); !done; curr, done = NextBed(file) {
		output <- curr
	}
	close(output)
}
