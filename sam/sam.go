package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"
)

type Sam struct {
	Header *SamHeader
	Aln    []*SamAln
}

type SamHeader struct {
	Text   []string
	Chroms []*chromInfo.ChromInfo
}

type SamAln struct {
	QName string
	Flag  int64
	RName string
	Pos   int64
	MapQ  int64 // mapping quality
	Cigar []*cigar.Cigar
	RNext string
	PNext int64
	TLen  int64
	Seq   []dna.Base
	Qual  string
	Extra string
}

func processHeaderLine(header *SamHeader, line string) {
	var err error
	var chrCount int64 = 0

	header.Text = append(header.Text, line)

	if strings.HasPrefix(line, "@SQ") && strings.Contains(line, "SN:") && strings.Contains(line, "LN:") {
		curr := chromInfo.ChromInfo{Name: "", Size: 0, Order: chrCount}
		chrCount++
		words := strings.Fields(line)
		for i := 1; i < len(words); i++ {
			elements := strings.Split(words[i], ":")
			switch elements[0] {
			case "SN":
				curr.Name = elements[1]
			case "LN":
				curr.Size, err = strconv.ParseInt(elements[1], 10, 64)
				if err != nil {
					log.Fatal(fmt.Errorf("Error: was expecting an integer length, but got:%s\n", elements[1]))
				}
			}
		}
		if curr.Name == "" || curr.Size == 0 {
			log.Fatal(fmt.Errorf("Thought I would get a name and non-zero size on this line: %s\n", line))
		}
		header.Chroms = append(header.Chroms, &curr)
	}
}

func ReadHeader(er *fileio.EasyReader) *SamHeader {
	var line string
	var err error
	var nextBytes []byte
	var header SamHeader

	for nextBytes, err = er.Peek(1); nextBytes[0] == '@' && err == nil; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		processHeaderLine(&header, line)
	}
	return &header
}

func processAlignmentLine(line string) *SamAln {
	var curr SamAln
	var err error

	words := strings.SplitN(line, "\t", 12)
	if len(words) < 11 {
		log.Fatal(fmt.Errorf("Was expecting atleast 11 columns per line, but this line did not:%s\n", line))
	}
	curr.QName = words[0]
	curr.Flag, err = strconv.ParseInt(words[1], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.RName = words[2]
	curr.Pos, err = strconv.ParseInt(words[3], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.MapQ, err = strconv.ParseInt(words[4], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.Cigar = cigar.FromString(words[5])
	curr.RNext = words[6]
	curr.PNext, err = strconv.ParseInt(words[7], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.TLen, err = strconv.ParseInt(words[8], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.Seq = dna.StringToBases(words[9])
	curr.Qual = words[10]
	if len(words) > 11 {
		curr.Extra = words[11]
	}
	return &curr
}

func NextAlignment(reader *fileio.EasyReader) (*SamAln, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	return processAlignmentLine(line), false
}

func ReadAlignments(er *fileio.EasyReader) []*SamAln {
	var line string
	var done bool
	var answer []*SamAln
	for line, done = fileio.EasyNextLine(er); !done; line, done = fileio.EasyNextLine(er) {
		answer = append(answer, processAlignmentLine(line))
	}
	return answer
}

func Read(filename string) (*Sam, error) {
	file := fileio.EasyOpen(filename)
	defer file.Close()

	header := ReadHeader(file)
	alnRecords := ReadAlignments(file)
	return &Sam{Header: header, Aln: alnRecords}, nil
}

func WriteHeaderToFileHandle(file *os.File, header *SamHeader) error {
	var err error

	for i, _ := range header.Text {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[i])
		common.ExitIfError(err)
	}
	return nil
}

func SamAlnToString(aln *SamAln) string {
	var answer string
	if aln.Extra == "" {
		answer = fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, cigar.ToString(aln.Cigar), aln.RNext, aln.PNext, aln.TLen, dna.BasesToString(aln.Seq), aln.Qual)
	} else {
		answer = fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s", aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, cigar.ToString(aln.Cigar), aln.RNext, aln.PNext, aln.TLen, dna.BasesToString(aln.Seq), aln.Qual, aln.Extra)
	}
	return answer
}

func WriteAlnToFileHandle(file *os.File, aln *SamAln) {
	_, err := fmt.Fprintf(file, "%s\n", SamAlnToString(aln))
	common.ExitIfError(err)
}

func SamChanToFile(incomingSams <-chan *SamAln, file *os.File, wg *sync.WaitGroup) {
	for alignedRead := range incomingSams {
		WriteAlnToFileHandle(file, alignedRead)
	}
	wg.Done()
}

func TestSamChanToFile(incomingSams <-chan *SamAln, file *os.File, wg *sync.WaitGroup) {

	for alignedRead := range incomingSams {
		WriteAlnToFileHandle(file, alignedRead)
	}
	wg.Done()
}
func SamChanToStdOut(incomingSams <-chan *SamAln, wg *sync.WaitGroup) {
	for alignedRead := range incomingSams {
		log.Printf("%s\n", SamAlnToString(alignedRead))
		//WriteAlnToFileHandle(file, alignedRead)

	}
	wg.Done()
}

func Write(filename string, data *Sam) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	err = WriteHeaderToFileHandle(file, data.Header)
	for i, _ := range data.Aln {
		WriteAlnToFileHandle(file, data.Aln[i])
	}
	return err
}

func AlignmentHeader(ref []*fasta.Fasta) *SamHeader {
	var header SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string
	for i := 0; i < len(ref); i++ {
		words = "@SQ\tSN:" + ref[i].Name + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})
	}
	return &header
}
