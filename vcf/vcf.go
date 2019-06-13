package vcf

import (
	"os"
	"io"
	"fmt"
	"bufio"
	"time"
	"strconv"
	"strings"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/chromInfo"
)

type Vcf struct {
	Chr     string
	Pos     int64
	Id      string
	Ref     string
	Alt     string
	Qual    float64
	Filter  string
	Info    string
	Format  string
	Unknown string
}

func Read(filename string) []*Vcf {
	var answer []*Vcf
	var curr *Vcf
	var line string
	file, err := os.Open(filename)
	if err != nil {
		return nil
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil
	}
	var err2 error
	var rline []byte
	for ; err2 != io.EOF; rline, _, err2 = reader.ReadLine() {
		line = string(rline[:])
		data := strings.Split(line, "\t")
		//fmt.Println("there is data here")
		switch {
		case strings.HasPrefix(line, "#"):
			//don't do anything
		case len(data) == 1:
			//these lines are sequences, and we are not recording them
			//fmt.Println("found sequences")
		case len(line) == 0:
			//blank line
			//fmt.Println("found blank")
		case len(data) == 10:
			//fmt.Println("found header line")
			position, _ := strconv.ParseInt(data[1], 10, 64)
			//qualFloat, _ := strconv.ParseFloat(data[5], 64)
			curr = &Vcf{Chr: data[0], Pos: position, Id: data[2], Ref: data[3], Alt: data[4], Qual: 0, Filter: data[6], Info: data[7], Format: data[8], Unknown: data[9]}
			answer = append(answer, curr)
		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer
}

//split vcf into slices to deal with different chromosomes
func VcfSplit(vcfRecord []*Vcf, fastaRecord []*fasta.Fasta) [][]*Vcf {
	var answer [][]*Vcf

	for i := 0; i < len(fastaRecord); i++ {
		var curr []*Vcf
		for j :=0; j < len(vcfRecord); j++ {
			var pointer *Vcf
			if strings.Compare(fastaRecord[i].Name, vcfRecord[j].Chr ) == 0 {
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
	header := BasicHeader()
	for h := 0; h < len(header); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header[h])
	}
	for i := 0; i < len(input); i++ {
		_, err = fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format, input[i].Unknown)
		common.ExitIfError(err)
	}
	return err
}

func Write(filename string, data []*Vcf) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	WriteVcfToFileHandle(file, data)
}

func PrintVcf(input []*Vcf) {
	for i := range input {
		fmt.Printf("%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format, input[i].Unknown)
	}
}

func PrintHeader(header []string) {
	for i := range header {
		fmt.Println(header[i])
	}
}

func ChromSizeHeader(chrom []*chromInfo.ChromInfo) []string{
	var header []string
	var line string
	t := time.Now()
	header = append(header, "##fileformat=VCFv4.2\n" +
		"##fileDate=" + t.Format("20060102") + "\n" +
		"##source=github.com/vertgenlab/gonomics\n" +
		"##reference=gasAcu1")
	for i := 0; i < len(chrom); i++ {
		line = "##contig=<ID=" + chrom[i].Name + ",length=" + string(chrom[i].Size) + ">"
		header = append(header, line)
	}
	header = append(header, "##phasing=none\n" +
		"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">\n" + 
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">" +
		"##INFO=<ID=QUERY,Number=1,Type=String,Description=\"Name of read/contig aligned\">" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n" +
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observation for each allele\">\n" +
		"##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n" +
		"##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">\n" +
		"##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n" +
		"##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">\n" +
		"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">\n" +
		"#CHROM  POS     ID      REF     ALT     QUAL    FILTER    INFO    FORMAT    Unknown")
	return header
}

func BasicHeader() []string{
	var header []string
	t := time.Now()
	header = append(header, "##fileformat=VCFv4.2\n" +
		"##fileDate=" + t.Format("20060102") + "\n" +
		"##source=github.com/vertgenlab/gonomics\n" +
		"##reference=gasAcu1")
	header = append(header, "##phasing=none\n" +
		"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">\n" + 
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">" +
		"##INFO=<ID=QUERY,Number=1,Type=String,Description=\"Name of read/contig aligned\">" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n" +
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observation for each allele\">\n" +
		"##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n" +
		"##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">\n" +
		"##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n" +
		"##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">\n" +
		"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">\n" +
		"#CHROM  POS     ID      REF     ALT     QUAL    FILTER    INFO    FORMAT    Unknown")

	return header
}