package sam

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type Sam struct {
	Header *SamHeader
	Aln    []*SamAln
}

type ChromSize struct {
	Name string
	Size int64
}

type SamHeader struct {
	Text      []string
	ChromInfo []*ChromSize
}

type SamAln struct {
	QName string
	Flag  int64
	RName string
	Pos   int64
	MapQ  int64  // mapping quality
	Cigar string // will change to align.Cigar
	RNext string
	PNext int64
	TLen  int64
	Seq   string
	Qual  string
	Extra string
}

/*
func samToVcf(alignment[]*Sam) {
	var answer []*vcf.Vcf
	for i := range alignment {

	}
}*/

func processHeaderLine(header *SamHeader, line string) {
	var err error

	header.Text = append(header.Text, line)

	if strings.HasPrefix(line, "@SQ") && strings.Contains(line, "SN:") && strings.Contains(line, "LN:") {
		curr := ChromSize{Name: "", Size: 0}
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
		header.ChromInfo = append(header.ChromInfo, &curr)
	}
}

func ReadHeader(reader *bufio.Reader) *SamHeader {
	var line string
	var err error
	var nextBytes []byte
	var header SamHeader

	for nextBytes, err = reader.Peek(1); nextBytes[0] == '@' && err == nil; nextBytes, err = reader.Peek(1) {
		line, err = reader.ReadString('\n')
		if err != nil {
			log.Fatal(err)
		}
		line = strings.TrimSuffix(line, "\n")
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
	curr.Cigar = words[5]
	curr.RNext = words[6]
	curr.PNext, err = strconv.ParseInt(words[7], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.TLen, err = strconv.ParseInt(words[8], 10, 64)
	if err != nil {
		log.Fatal(err)
	}
	curr.Seq = words[9]
	curr.Qual = words[10]
	if len(words) > 11 {
		curr.Extra = words[11]
	}
	return &curr
}

func ReadAlignments(reader *bufio.Reader) []*SamAln {
	var line string
	var err error
	var answer []*SamAln
	for line, err = reader.ReadString('\n'); err == nil; line, err = reader.ReadString('\n') {
		line = strings.TrimSuffix(line, "\n")
		answer = append(answer, processAlignmentLine(line))
	}
	if err != io.EOF {
		log.Fatal(err)
	}
	return answer
}

func Read(filename string) (*Sam, error) {
	file, err := os.Open(filename)
	defer file.Close()
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)

	header := ReadHeader(reader)
	alnRecords := ReadAlignments(reader)
	return &Sam{Header: header, Aln: alnRecords}, nil
}

func WriteHeaderToFileHandle(file *os.File, header *SamHeader) error {
	var err error

	for i, _ := range header.Text {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[i])
		if err != nil {
			return err
		}
	}
	return nil
}

func WriteAlnToFileHandle(file *os.File, aln *SamAln) error {
	var err error
	if aln.Extra == "" {
		_, err = fmt.Fprintf(file, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n", aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, aln.Cigar, aln.RNext, aln.PNext, aln.TLen, aln.Seq, aln.Qual)
	} else {
		_, err = fmt.Fprintf(file, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n", aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, aln.Cigar, aln.RNext, aln.PNext, aln.TLen, aln.Seq, aln.Qual, aln.Extra)
	}
	return err
}

func Write(filename string, data *Sam) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	err = WriteHeaderToFileHandle(file, data.Header)
	for i, _ := range data.Aln {
		err = WriteAlnToFileHandle(file, data.Aln[i])
		if err != nil {
			return err
		}
	}
	return err
}

/*
func Read(filename string) ([]*Sam, error) {
	var answer []*Sam
	var curr *Sam
	var line string
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil, err
	}
	var err2 error
	var rline []byte
	for ;err2 != io.EOF; rline, _, err2 = reader.ReadLine()  {
		line = string(rline[:])

		if strings.HasPrefix(line, "@") {
			fmt.Println(line)
		}
		data := strings.Split(line, "\t")

		switch {
		case len(data) == 10:
			bf, _ := strconv.ParseInt(data[1], 10, 64)
			qf, _ := strconv.ParseInt(data[4], 10, 64)
			pn, _ := strconv.ParseInt(data[7], 10, 64)
			tl, _ := strconv.ParseInt(data[8], 10, 64)
			//curr = &Sam{Chr: data[0], Pos: position, Id: data[2], Ref: data[3], Alt: data[4], Qual: 0, Filter: data[6], Info: data[7], Format: data[8], Unknown: data[9]}
			curr = &Sam{QName: data[0], BitFlag: bf, RefName: data[2], CurrPos: data[3], QMap: qf, Cigar: data[5], RNext: data[6], PosNext: pn, TmpLen: tl, Seq: data[9], QualBase: " "}
			answer = append(answer, curr)

		case len(data) == 11:
			bf, _ := strconv.ParseInt(data[1], 10, 64)
			qf, _ := strconv.ParseInt(data[4], 10, 64)
			pn, _ := strconv.ParseInt(data[7], 10, 64)
			tl, _ := strconv.ParseInt(data[8], 10, 64)
			//curr = &Sam{Chr: data[0], Pos: position, Id: data[2], Ref: data[3], Alt: data[4], Qual: 0, Filter: data[6], Info: data[7], Format: data[8], Unknown: data[9]}
			curr = &Sam{QName: data[0], BitFlag: bf, RefName: data[2], CurrPos: data[3], QMap: qf, Cigar: data[5], RNext: data[6], PosNext: pn, TmpLen: tl, Seq: data[9], QualBase: data[10]}
			answer = append(answer, curr)
		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer, nil
}



func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	samFile := flag.Arg(0)
	//out, _ := ReadIn(samFile)
	align, _ := ReadIn(samFile)

	for i := range align {
			fmt.Println(align[i])

	}
	//fmt.Println(len(align))
}
*/
//linebyline := metaData{qName: colData[0], bitFlag: bf, refName: colData[2], currPos: colData[3], qMap: qf, cigar: colData[5], rNext: colData[6], posNext: pn, tmpLen: tl, seq: colData[9], qualBase: colData[10]}
//currSam := Sam{Header: tmpHeader, Alignment: &linebyline}
