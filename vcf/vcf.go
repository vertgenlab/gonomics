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
	"time"
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

//Might get rid of this
type VCF struct {
	Header *VcfHeader
	Vcf    []*Vcf
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

func GoReadToChan(filename string) (<-chan *Vcf, *VcfHeader) {
	output := make(chan *Vcf)
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	go ReadToChan(file, output)
	return output, header
}

func ReadToChan(file *fileio.EasyReader, output chan<- *Vcf) {
	var curr *Vcf
	var done bool
	for curr, done = NextVcf(file); !done; curr, done = NextVcf(file) {
		output <- curr
	}
	close(output)
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

func ReadVcf(filename string) *VCF {
	file := fileio.EasyOpen(filename)
	defer file.Close()

	header := ReadHeader(file)
	vcfRecords := Read(filename)
	return &VCF{Header: header, Vcf: vcfRecords}
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

func processHeader(header *VcfHeader, line string) {
	var err error
	if strings.HasPrefix(line, "#") {
		header.Text = append(header.Text, line)
	}
	if err != nil {
		log.Fatal("There was an error reading the header line")
	}
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

func WriteHeader(file *os.File) {
	var err error
	header := MakeHeader()
	for h := 0; h < len(header); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header[h])
	}
	common.ExitIfError(err)
}

func WriteBetterHeader(file *fileio.EasyWriter, fa []*fasta.Fasta) {
	var err error
	header := BetterHeader(fa)
	for h := 0; h < len(header); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header[h])
	}
	common.ExitIfError(err)
}

func NewWrite(filename string, data []*Vcf, fa []*fasta.Fasta) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteBetterHeader(file, fa)
	for _, each := range data {
		WriteVcf(file, each)
	}
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
	WriteHeader(file)
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
	var err error
	file := fileio.EasyCreate("/dev/stdout")
	defer file.Close()
	_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", data.Chr, data.Pos, data.Id, data.Ref, data.Alt, data.Qual, data.Filter, data.Info, data.Format)
	common.ExitIfError(err)
}

func PrintHeader(header []string) {
	for i := range header {
		fmt.Println(header[i])
	}
}

//TODO will be removed in the next few weeks
func MakeHeader() []string {
	//TODO: add logic to add contig length to header of file
	var header []string
	//var line string
	t := time.Now()
	header = append(header, "##fileformat=VCFv4.2\n"+
		"##fileDate="+t.Format("20060102")+"\n"+
		"##source=github.com/vertgenlab/gonomics")
	header = append(header, "##phasing=none\n"+
		"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">\n"+
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">"+
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">"+
		"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"+
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"+
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"+
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observation for each allele\">\n"+
		"##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n"+
		"##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">\n"+
		"##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n"+
		"##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">\n"+
		"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">\n"+
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNOTES")
	return header
}

func BetterHeader(fa []*fasta.Fasta) []string {
	var header []string
	t := time.Now()
	header = append(header, "##fileformat=VCFv4.2\n"+
		"##fileDate="+t.Format("20060102")+"\n"+
		"##source=github.com/vertgenlab/gonomics")
	header = append(header, "##phasing=none\n"+
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">"+
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">"+
		"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"+
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	var text string = ""
	if fa != nil {
		for i := 0; i < len(fa); i++ {
			text += fmt.Sprintf("##contig=<ID=%s,length=%d>\n", fa[i].Name, len(fa[i].Seq))
		}
	}
	text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNOTES"
	header = append(header, text)
	return header
}
