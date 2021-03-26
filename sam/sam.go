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
	"strconv"
	"strings"
	"sync"
)

type Sam struct {
	Header Header
	Aln    []Aln
}

type Aln struct {
	QName string         // query template name: [!-?A-~]{1,254}
	Flag  uint16         // bitwise flag, bits defined in info.go
	MapQ  uint8          // mapping quality
	RName string         // reference sequence name: \*|[:rname:∧ *=][:rname:]*
	Pos   int            // 1-based leftmost mapping position
	Cigar []*cigar.Cigar // parsed cigar: originally string with \*|([0-9]+[MIDNSHPX=])+
	RNext string         // reference name of the mate/next read: \*|=|[:rname:∧ *=][:rname:]*
	PNext int            // position of the mate/next read
	TLen  int            // observed template length
	Seq   []dna.Base     // parsed sequence: originally string with \*|[A-Za-z=.]+
	Qual  string         // ASCII of Phred-scaled base quality+33: [!-~]+
	Extra string         // TODO parse to map or slice w/ index embedded in header???
}

func ReadToChan(file *fileio.EasyReader, data chan<- Aln, wg *sync.WaitGroup) {
	for curr, done := NextAlignment(file); !done; curr, done = NextAlignment(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

func GoReadToChan(filename string) (<-chan Aln, Header) {
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	var wg sync.WaitGroup
	data := make(chan Aln)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, header
}

func SamChanToFile(incomingSams <-chan Aln, filename string, header Header, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if header.Text != nil {
		WriteHeaderToFileHandle(file, header)
	}
	for alignedRead := range incomingSams {
		WriteAlnToFileHandle(file, alignedRead)
	}
	wg.Done()
}

func processAlignmentLine(line string) Aln {
	var curr Aln
	var err error
	var currUint uint64

	words := strings.SplitN(line, "\t", 12)
	if len(words) < 11 {
		log.Fatal(fmt.Errorf("Was expecting atleast 11 columns per line, but this line did not:%s\n", line))
	}
	curr.QName = words[0]
	currUint, err = strconv.ParseUint(words[1], 10, 16)
	if err != nil {
		log.Fatal(err)
	}
	curr.Flag = uint16(currUint)
	curr.RName = words[2]
	curr.Pos, err = strconv.Atoi(words[3])
	if err != nil {
		log.Fatal(err)
	}
	currUint, err = strconv.ParseUint(words[4], 10, 8)
	if err != nil {
		log.Fatal(err)
	}
	curr.MapQ = uint8(currUint)
	curr.Cigar = cigar.FromString(words[5])
	curr.RNext = words[6]
	curr.PNext, err = strconv.Atoi(words[7])
	if err != nil {
		log.Fatal(err)
	}
	curr.TLen, err = strconv.Atoi(words[8])
	if err != nil {
		log.Fatal(err)
	}
	curr.Seq = dna.StringToBases(words[9])
	curr.Qual = words[10]
	if len(words) > 11 {
		curr.Extra = words[11]
	}
	return curr
}

func NextAlignment(reader *fileio.EasyReader) (Aln, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return Aln{}, true
	}
	return processAlignmentLine(line), false
}

func ReadAlignments(er *fileio.EasyReader) []Aln {
	var line string
	var done bool
	var answer []Aln
	for line, done = fileio.EasyNextLine(er); !done; line, done = fileio.EasyNextLine(er) {
		answer = append(answer, processAlignmentLine(line))
	}
	return answer
}

func Read(filename string) (Sam, error) {
	file := fileio.EasyOpen(filename)
	defer file.Close()

	header := ReadHeader(file)
	alnRecords := ReadAlignments(file)
	return Sam{Header: header, Aln: alnRecords}, nil
}

func SamAlnToString(aln Aln) string {
	var answer string
	if aln.Extra == "" {
		answer = fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, cigar.ToString(aln.Cigar), aln.RNext, aln.PNext, aln.TLen, dna.BasesToString(aln.Seq), aln.Qual)
	} else {
		answer = fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s", aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, cigar.ToString(aln.Cigar), aln.RNext, aln.PNext, aln.TLen, dna.BasesToString(aln.Seq), aln.Qual, aln.Extra)
	}
	return answer
}

func WriteAlnToFileHandle(file *fileio.EasyWriter, aln Aln) {
	_, err := fmt.Fprintf(file, "%s\n", SamAlnToString(aln))
	common.ExitIfError(err)
}

func Write(filename string, data Sam) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteHeaderToFileHandle(file, data.Header)
	for i := range data.Aln {
		WriteAlnToFileHandle(file, data.Aln[i])
	}
}

func ReadHeader(er *fileio.EasyReader) Header {
	var line string
	var err error
	var nextBytes []byte
	var header Header

	for nextBytes, err = er.Peek(1); err == nil && nextBytes[0] == '@'; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		processHeaderLine(header, line)
	}
	return header
}

func processHeaderLine(header Header, line string) {
	var chrCount int = 0

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
				curr.Size = common.StringToInt(elements[1])
			}
		}
		if curr.Name == "" || curr.Size == 0 {
			//	log.Fatal(fmt.Errorf("Thought I would get a name and non-zero size on this line: %s\n", line))
		}
		header.Chroms = append(header.Chroms, curr)
	}
}

func WriteHeaderToFileHandle(file *fileio.EasyWriter, header Header) error {
	var err error

	for i := range header.Text {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[i])
		common.ExitIfError(err)
	}
	return nil
}

func ChromInfoSamHeader(chromSize []chromInfo.ChromInfo) Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(chromSize); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", chromSize[i].Name, chromSize[i].Size)
		header.Text = append(header.Text, words)
	}
	return header
}

func ChromInfoMapSamHeader(chromSize map[string]chromInfo.ChromInfo) *Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string
	var i int
	for i = 0; i < len(chromSize); {
		for j := range chromSize {
			if i == chromSize[j].Order {
				words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", chromSize[j].Name, chromSize[j].Size)
				header.Text = append(header.Text, words)
				i++
			}
		}
	}
	return &header
}

func FastaHeader(ref []fasta.Fasta) Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", ref[i].Name, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, chromInfo.ChromInfo{Name: ref[i].Name, Size: len(ref[i].Seq)})
	}
	return header
}
