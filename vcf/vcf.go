package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"os"
	"strings"
	"sync"
)

//slices of slice
type Vcf struct {
	Chr    string
	Pos    int64
	Id     string
	Ref    string
	Alt    string
	Qual   float64
	Filter string
	Info   string
	Format string
	Notes  string
}

type VcfHeader struct {
	Text []string
}

func Read(filename string) []*Vcf {
	var answer []*Vcf
	file := fileio.EasyOpen(filename)
	defer file.Close()
	ReadHeader(file)
	var curr *Vcf
	var done bool
	for curr, done = NextVcf(file); !done; curr, done = NextVcf(file) {
		answer = append(answer, curr)
	}
	return answer
}

func ReadToChan(file *fileio.EasyReader, data chan<- *Vcf, wg *sync.WaitGroup) {
	for curr, done := NextVcf(file); !done; curr, done = NextVcf(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

func GoReadToChan(filename string) (<-chan *Vcf, *VcfHeader) {
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	var wg sync.WaitGroup
	data := make(chan *Vcf)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, header
}

func processVcfLine(line string) *Vcf {
	var curr *Vcf
	data := strings.SplitN(line, "\t", 10)
	//switch {
	//case strings.HasPrefix(line, "#"):
	//don't do anything
	//case len(data) == 1:
	//these lines are sequences, and we are not recording them
	//case len(line) == 0:
	//blank line
	if len(data) < 9 {
		log.Fatalf("Error when reading this vcf line:\n%s\nExpecting at least 9 columns", line)
	}
	curr = &Vcf{Chr: data[0], Pos: common.StringToInt64(data[1]), Id: data[2], Ref: data[3], Alt: data[4], Filter: data[6], Info: data[7], Format: data[8], Notes: ""}
	if strings.Compare(data[5], ".") == 0 {
		curr.Qual = 255
	} else {
		curr.Qual = common.StringToFloat64(data[5])
	}
	if len(data) > 9 {
		curr.Notes = data[9]
	}
	return curr
}

func NextVcf(reader *fileio.EasyReader) (*Vcf, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return nil, true
	}
	return processVcfLine(line), false
}

func ReadHeader(er *fileio.EasyReader) *VcfHeader {
	var line string
	var err error
	var nextBytes []byte
	var header VcfHeader
	for nextBytes, err = er.Peek(1); err == nil && nextBytes[0] == '#'; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		processHeader(&header, line)
	}
	return &header
}

//split vcf into slices to deal with different chromosomes
func VcfSplit(vcfRecord []*Vcf, fastaRecord []*fasta.Fasta) [][]*Vcf {
	var answer [][]*Vcf
	for i := 0; i < len(fastaRecord); i++ {
		var curr []*Vcf
		for j := 0; j < len(vcfRecord); j++ {
			var pointer *Vcf
			if strings.Compare(fastaRecord[i].Name, vcfRecord[j].Chr) == 0 {
				pointer = vcfRecord[j]
				curr = append(curr, pointer)
			}
		}
		Sort(curr)
		answer = append(answer, curr)
	}
	return answer
}

//TODO(craiglowe): Look into unifying WriteVcfToFileHandle and WriteVcf and benchmark speed
func WriteVcfToFileHandle(file *os.File, input []*Vcf) {
	var err error
	for i := 0; i < len(input); i++ {
		if input[i].Notes == "" {
			_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format)
		} else {
			_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format, input[i].Notes)
		}
		common.ExitIfError(err)
	}
}

func WriteVcf(file io.Writer, input *Vcf) {
	var err error
	if input.Notes == "" {
		_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", input.Chr, input.Pos, input.Id, input.Ref, input.Alt, input.Qual, input.Filter, input.Info, input.Format)
	} else {
		_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input.Chr, input.Pos, input.Id, input.Ref, input.Alt, input.Qual, input.Filter, input.Info, input.Format, input.Notes)
	}
	common.ExitIfError(err)
}

func Write(filename string, data []*Vcf) {
	file := fileio.MustCreate(filename)
	defer file.Close()

	WriteVcfToFileHandle(file, data)
}

func PrintVcf(data []*Vcf) {
	for i := 0; i < len(data); i++ {
		PrintSingleLine(data[i])
	}
}

func PrintVcfLines(data []*Vcf, num int) {
	for i := 0; i < num; i++ {
		PrintSingleLine(data[i])
	}
}

func PrintSingleLine(data *Vcf) {
	fmt.Printf("%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", data.Chr, data.Pos, data.Id, data.Ref, data.Alt, data.Qual, data.Filter, data.Info, data.Format)
}

//Checks suffix of filename to confirm if the file is a vcf formatted file
func IsVcfFile(filename string) bool {
	if strings.HasSuffix(filename, ".vcf") || strings.HasSuffix(filename, ".vcf.gz") {
		return true
	} else {
		return false
	}
}
