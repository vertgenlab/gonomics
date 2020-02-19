package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
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

type VCF struct {
	Header *VcfHeader
	Vcf    []*Vcf
}

type VcfHeader struct {
	Text []string
}

func processVcfLine(line string) *Vcf {

	data := Vcf{Chr: "", Pos: 0, Id: "", Ref: "", Alt: "", Qual: 0, Filter: "", Info: "", Format: "", Notes: ""}
	text := strings.Split(line, "\t")
	if len(text) < 9 {
		log.Fatal(fmt.Errorf("Was expecting atleast 8 columns per line, but this line did not:%s\n", line))
	}
	data.Pos = common.StringToInt64(text[1])
	data.Qual = common.StringToFloat64(text[5])

	data.Chr = text[0]
	data.Id = text[2]
	data.Ref = text[3]
	data.Alt = text[4]
	data.Filter = text[6]
	data.Info = text[7]
	data.Format = text[8]
	if len(text) > 9 {
		data.Notes = text[9]
	}

	return &data
}

func NextVcf(reader *fileio.EasyReader) (*Vcf, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	return processVcfLine(line), false
}

func Read(filename string) []*Vcf {
	file := fileio.EasyOpen(filename)
	defer file.Close()

	var line string
	var done bool
	var answer []*Vcf
	ReadHeader(file)
	for line, done = fileio.EasyNextLine(file); !done; line, done = fileio.EasyNextLine(file) {
		answer = append(answer, processVcfLine(line))
	}
	return answer
}

func ReadVcf(filename string) *VCF {
	file := fileio.EasyOpen(filename)
	defer file.Close()

	header := ReadHeader(file)
	vcfRecords := Read(filename)
	return &VCF{Header: header, Vcf: vcfRecords}
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

func ReadHeader(er *fileio.EasyReader) *VcfHeader {
	var line string
	var err error
	var nextBytes []byte
	var header VcfHeader

	for nextBytes, err = er.Peek(1); nextBytes[0] == '#' && err == nil; nextBytes, err = er.Peek(1) {
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

func WriteVcfToFileHandle(file *os.File, input []*Vcf) error {
	var err error

	//var chrLen []*chromInfo.ChromInfo = []*chromInfo.ChromInfo{}
	header := MakeHeader()
	for h := 0; h < len(header); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header[h])
	}
	for i := 0; i < len(input); i++ {
		if input[i].Notes == "" {
			_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format)

		} else {
			_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format, input[i].Notes)
		}

	}
	common.ExitIfError(err)
	return err
}

func Write(filename string, data []*Vcf) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	WriteVcfToFileHandle(file, data)
}

func PrintVcf(data []*Vcf) {
	Write("/dev/stdout", data)
}

func PrintHeader(header []string) {
	for i := range header {
		fmt.Println(header[i])
	}
}

func MakeHeader() []string {
	//TODO: add logic to add contig length to header of file
	var header []string
	//var line string
	t := time.Now()
	header = append(header, "##fileformat=VCFv4.2\n"+
		"##fileDate="+t.Format("20060102")+"\n"+
		"##source=github.com/vertgenlab/gonomics\n"+
		"##reference=gasAcu1")
	//if len(chrom) > 0 {
	//	for i := 0; i < len(chrom); i++ {
	//		line = "##contig=<ID=" + chrom[i].Name + ",length=" + string(chrom[i].Size) + ">"
	//		header = append(header, line)
	//}
	//}
	header = append(header, "##phasing=none\n"+
		"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">\n"+
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">"+
		"##INFO=<ID=QUERY,Number=1,Type=String,Description=\"Name of read/contig aligned\">"+
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
