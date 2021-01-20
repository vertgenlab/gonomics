//Package vcf contains functions for reading, writing, and manipulating VCF format files. More information on the VCF file format can be found
//in its official documentation at https://samtools.github.io/hts-specs/VCFv4.2.pdf. This file is parsed into a struct containing header information
//as well as a Vcf struct containing the information from each data line.

package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
	"strconv"
	"sync"
)

//Vcf contains information for each line of a VCF format file, corresponding to variants at one position of a reference genome.
type Vcf struct {
	Chr     string
	Pos     int
	Id      string
	Ref     string
	Alt     []string //TODO: bytes.buffer for memory optimization
	Qual    float64
	Filter  string
	Info    string
	Format  []string
	Samples []GenomeSample
}

//GenomeSample is a substruct of Vcf, and contains information about each individual represented in a VCF line.
// AlleleOne and AlleleTwo are set to -1 if no genotype data is present.
type GenomeSample struct {
	AlleleOne  int16    //First allele in genotype, 0 for reference, 1 for Alt[0], 2 for Alt[1], etc.
	AlleleTwo  int16    //Second allele in genotype, same number format as above.
	Phased     bool     //True for phased genotype, false for unphased.
	FormatData []string //FormatData contains additional sample fields after the genotype, which are parsed into a slice delimited by colons. Currently contains a dummy empty strin in FormatData[0] corresponding to "GT" in Format, so indices in FormatData will match the indices in Format.
}

//VcfHeader contains all of the information present in the header section of a VCf, simply delimited by line to form a slice of strings.
type VcfHeader struct {
	Text []string
}

//Read parses a slice of VCF structs from an input filename. Does not store or return the header.
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

//ReadToChan is a helper function of GoReadToChan.
func ReadToChan(file *fileio.EasyReader, data chan<- *Vcf, wg *sync.WaitGroup) {
	for curr, done := NextVcf(file); !done; curr, done = NextVcf(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

//GoReadToChan parses VCF structs from an input filename and returns a chan of VCF structs along with the header of the VCF file.
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

//processVcfLine is a helper function of NextVcf that parses a VCF struct from an input line of a VCF file provided as a string.
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
	curr = &Vcf{Chr: data[0], Pos: common.StringToInt(data[1]), Id: data[2], Ref: data[3], Alt: strings.Split(data[4], ","), Filter: data[6], Info: data[7], Format: strings.Split(data[8], ":")}
	if strings.Compare(data[5], ".") == 0 {
		curr.Qual = 255
	} else {
		curr.Qual = common.StringToFloat64(data[5])
	}
	if len(data) > 9 {
		//DEBUG: fmt.Println(line)
		curr.Samples = ParseNotes(data[9], curr.Format)
	}
	return curr
}

//ParseNotes is a helper function of processVcfLine. Generates a slice of GenomeSample structs from a VCF data line.
func ParseNotes(data string, format []string) []GenomeSample {
	//DEBUG: fmt.Printf("Format: %s. Format[0]: %s.\n", format, format[0])
	if len(format) == 0 {
		log.Fatalf("Parsing error: cannot parse notes without formatting information.")
	}
	//firstFormat := format[0]
	//fmt.Println(firstFormat)
	if format[0] != "GT" && format[0] != "." { //len(format) == 0 checks for when there is info in the notes column but no format specification
		log.Fatalf("VCF format files with sample information must begin with \"GT\" as the first format column or be marked blank with a period. Here was the first format entry: %s.\nError parsing the line with this Notes column: %s.\n", format[0], data)
	}

	if format[0] == "." {//if the format column is blank, we do not need to parse further.
		var blankAnswer []GenomeSample = make([]GenomeSample, 0)
		return blankAnswer
	}

	text := strings.Split(data, "\t")
	var fields []string
	var alleles []string
	var err error
	var n int64
	var answer []GenomeSample = make([]GenomeSample, len(text))
	for i := 0; i < len(text); i++ {
		fields = strings.Split(text[i], ":")
		if strings.Compare(fields[0], "./.") == 0 || strings.Compare(fields[0], ".|.") == 0 {
			answer[i] = GenomeSample{AlleleOne: -1, AlleleTwo: -1, Phased: false, FormatData: fields}
		} else if strings.Contains(fields[0], "|") {
			alleles = strings.SplitN(fields[0], "|", 2)
			answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: true, FormatData: fields}
		} else if strings.Contains(fields[0], "/") {
			alleles = strings.SplitN(fields[0], "/", 2)
			answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: false, FormatData: fields}
		} else {
			n, err = strconv.ParseInt(fields[0], 10, 16)
			if err == nil && n < int64(len(text)) {
				answer[i] = GenomeSample{AlleleOne: int16(n), AlleleTwo: -1, Phased: false, FormatData: fields}
			} else {
				log.Fatalf("Error: Unexpected parsing error on the following line:\n%s", data)
			}
		}
		answer[i].FormatData[0] = "" //clears the genotype from the first other slot, making this a dummy position
	}
	return answer
}

//NextVcf is a helper function of Read and GoReadToChan. Checks a reader for additional data lines and parses a Vcf line if more lines exist.
func NextVcf(reader *fileio.EasyReader) (*Vcf, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return nil, true
	}
	return processVcfLine(line), false
}

//ReadHeader is a helper function of GoReadToChan. Parses a VCF header with a reader.
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

//FormatToString converts the []string Format struct into a string by concatenating with a colon delimiter.
func FormatToString(format []string) string {
	if len(format) == 0 {
		return ""
	}
	var answer = format[0]
	for i := 1; i < len(format); i++ {
		answer = answer + ":" + format[i]
	}
	return answer
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

//TODO(craiglowe): Look into unifying WriteVcfToFileHandle and WriteVcf and benchmark speed. geno bool variable determines whether to print notes or genotypes.
func WriteVcfToFileHandle(file io.Writer, input []*Vcf) {
	var err error
	for i := 0; i < len(input); i++ {
		//DEBUG:fmt.Printf("Notes: %s\n", input[i].Notes)
		if len(input[i].Format) == 0 {
			_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, strings.Join(input[i].Alt, ","), input[i].Qual, input[i].Filter, input[i].Info)
		} else {
			_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, strings.Join(input[i].Alt, ","), input[i].Qual, input[i].Filter, input[i].Info, FormatToString(input[i].Format), SamplesToString(input[i].Samples))
		}
		common.ExitIfError(err)
	}
}

//WriteVcf writes an individual Vcf struct to an io.Writer.
func WriteVcf(file io.Writer, input *Vcf) {
	var err error
	if len(input.Format) == 0 {
		_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\n", input.Chr, input.Pos, input.Id, input.Ref, strings.Join(input.Alt, ","), input.Qual, input.Filter, input.Info)
	} else {
		_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input.Chr, input.Pos, input.Id, input.Ref, strings.Join(input.Alt, ","), input.Qual, input.Filter, input.Info, FormatToString(input.Format), SamplesToString(input.Samples))
	}
	common.ExitIfError(err)
}

//Write writes a []*Vcf to an output filename.
func Write(filename string, data []*Vcf) {
	file := fileio.MustCreate(filename)
	defer file.Close()

	WriteVcfToFileHandle(file, data)
}

//PrintVcf prints every line of a []*Vcf.
func PrintVcf(data []*Vcf) {
	for i := 0; i < len(data); i++ {
		PrintSingleLine(data[i])
	}
}

//PrintVcfLines prints the first "num" lines from a []*Vcf.
func PrintVcfLines(data []*Vcf, num int) {
	for i := 0; i < num; i++ {
		PrintSingleLine(data[i])
	}
}

//PrintSingleLine prints an individual Vcf line.
func PrintSingleLine(data *Vcf) {
	fmt.Printf("%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", data.Chr, data.Pos, data.Id, data.Ref, strings.Join(data.Alt, ","), data.Qual, data.Filter, data.Info, data.Format, SamplesToString(data.Samples))
}

//IsVcfFile checks suffix of filename to confirm if the file is a vcf formatted file
func IsVcfFile(filename string) bool {
	if strings.HasSuffix(filename, ".vcf") || strings.HasSuffix(filename, ".vcf.gz") {
		return true
	} else {
		return false
	}
}
