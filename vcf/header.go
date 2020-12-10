package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"io"
	"log"
	"strings"
	"time"
)

func processHeader(header *VcfHeader, line string) {
	if strings.HasPrefix(line, "#") {
		header.Text = append(header.Text, line)
	} else {
		log.Fatal("There was an error reading the header line")
	}
}

//If you have multiple samples to add to the header, use strings.Join(samples[], "\t") as an argument that combines multiple samples by tabs
//Name input is strictly used to push in a name for the sample column.

func NewHeader(name string) *VcfHeader {
	var header *VcfHeader = &VcfHeader{}
	t := time.Now()
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, "##fileDate="+t.Format("20060102"))
	header.Text = append(header.Text, "##source=github.com/vertgenlab/gonomics")
	header.Text = append(header.Text, "##phasing=none")
	header.Text = append(header.Text, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">")
	header.Text = append(header.Text, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">")
	header.Text = append(header.Text, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">")
	header.Text = append(header.Text, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", name))
	return header
}

//The other option is add everything as one string, so we don';t have to keep appending, parsing will in general be the same, but append is more consistant with how we read in line, by line.
/*
var text string =
		"##fileformat=VCFv4.2\n"+
		"##fileDate="+t.Format("20060102")+"\n"+
		"##source=github.com/vertgenlab/gonomics\n"+
		"##phasing=none\n"+
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">\n"+
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n"+
		"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"+
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"+
		fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", name)

*/

func AddContigHeader(header *VcfHeader, fa []*fasta.Fasta) *VcfHeader {
	var text string = ""
	for i := 0; i < len(fa); i++ {
		text += fmt.Sprintf("##contig=<ID=%s,length=%d>\n", fa[i].Name, len(fa[i].Seq))
	}
	header.Text = append(header.Text, text)
	return header
}

func NewWriteHeader(file io.Writer, header *VcfHeader) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
		common.ExitIfError(err)
	}
}

func WriteMultiSamplesHeader(file io.Writer, header *VcfHeader, listNames []string) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		if strings.Contains(header.Text[h], "#CHROM\t") {
			name := strings.Join(listNames, "\t")
			_, err = fmt.Fprintf(file, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", name))
			common.ExitIfError(err)
		} else {
			_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
			common.ExitIfError(err)
		}
	}
}

func PrintHeader(header *VcfHeader) {
	for i := 0; i < len(header.Text); i++ {
		fmt.Println(header.Text[i])
	}
}
