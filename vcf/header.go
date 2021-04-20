package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"io"
	"log"
	"strings"
	"time"
)

func processHeader(header Header, line string) Header {
	if strings.HasPrefix(line, "#") {
		header.Text = append(header.Text, line)
	} else {
		log.Fatal("There was an error reading the header line")
	}
	return header
}

//If you have multiple samples to add to the header, use strings.Join(samples[], "\t") as an argument that combines multiple samples by tabs
//Name input is strictly used to push in a name for the sample column.

func NewHeader(name string) Header {
	var header Header
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

func NewWriteHeader(file io.Writer, header Header) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
		common.ExitIfError(err)
	}
}

func WriteMultiSamplesHeader(file io.Writer, header Header, listNames []string) {
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

//Uses Vcf header to create 2 hash maps 1) is the sample index that maps the which allele each sample has in Vcf 2) hash reference chromsome names to an index (used to build uint64 containing chromID and position)
func HeaderToMaps(header Header) *SampleHash {
	var name string
	var index, hapIdx int16
	var hash *SampleHash = &SampleHash{Fa: make(map[string]int16), GIndex: make(map[string]int16)}
	for _, line := range header.Text {
		if strings.HasPrefix(line, "##contig") {
			name = strings.Split(strings.Split(line, "=")[2], ",")[0]
			_, ok := hash.Fa[name]
			if !ok {
				hash.Fa[name] = index
				index++
			}
		} else if strings.HasPrefix(line, "#CHROM") {
			words := strings.Split(line, "\t")[9:]
			for hapIdx = 0; hapIdx < int16(len(words)); hapIdx++ {
				hash.GIndex[words[hapIdx]] = hapIdx
			}
		}
	}
	return hash
}

//HeaderGetSampleList returns an ordered list of the samples present in the header of a Vcf file. Useful when adding or removing samples from a VCF.
func HeaderGetSampleList(header Header) []string {
	var answer []string
	for _, line := range header.Text {
		if strings.HasPrefix(line, "#CHROM") {
			return strings.Split(line, "\t")[9:]
		}
	}
	log.Fatalf("No Sample info in VCF line, cannot parse sample names.")
	return answer
}

//HeaderUpdateSampleList can be provided with a new list of samples to update the sample list in a Header.
func HeaderUpdateSampleList(header Header, newSamples []string) {
	var line string
	for i := 0; i < len(header.Text); i++ {
		if strings.HasPrefix(header.Text[i], "#CHROM") {
			line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
			for j := 0; j < len(newSamples); j++ {
				line += "\t" + newSamples[j]
			}
			header.Text[i] = line
		}
	}
}
