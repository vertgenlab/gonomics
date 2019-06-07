package vcf

import (
	"bufio"

	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"strconv"
	"strings"
	"io"
	"sort"
	"time"
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

type Snp struct {
	RefName   string
	RefPos    int64
	RefSub    string
	QueryName string
	QueryPos  int64
	QuerySub  string
	Next *Snp
}

func Read(filename string) []*Vcf {
	var answer []*Vcf
	var curr *Vcf
	var line string
	file, _ := os.Open(filename)

	reader := bufio.NewReader(file)
	
	var err2 error
	var rline []byte
	for ; err2 != io.EOF; rline, _, err2 = reader.ReadLine() {
		line = string(rline[:])
		data := strings.Split(line, "\t")
		//fmt.Println("there is data here")

		switch {
		case strings.HasPrefix(line, "#"):
			//don't do anything
			//fmt.Println("found #")
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
			curr = &Vcf{Chr: data[0], Pos: position, Id: data[2], Ref: data[3], Alt: data[4], Qual: 0, Filter: data[6], Info: data[7], Format: data[8]}
			answer = append(answer, curr)
		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer
}

func SnpsToVcf(snps []*Snp, reference map[string][]dna.Base, contigs map[string][]dna.Base) []*Vcf {
	var answer []*Vcf
	var curr *Vcf
	var currentSeq []dna.Base
	var contigPointer []dna.Base 
	for i := 0; i < len(snps); i++ {
		var infoTag string
		switch {
		//SNP
		case SnpTruth(snps[i]) == 0:
			infoTag = "POS=" + strconv.FormatInt(snps[i].QueryPos, 10)
			curr = &Vcf{Chr: snps[i].RefName, Pos: snps[i].RefPos, Id: snps[i].QueryName, Ref: snps[i].RefSub, Alt: snps[i].QuerySub, Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=SNP"}
			//fmt.Println(snps[i].RefSub, snps[i].QuerySub)
			answer = append(answer, curr)
		//logic for insertion in vcf record
		case strings.Compare(snps[i].RefSub, "-") == 0:
			var altTmp string
			currentSeq = reference[snps[i].RefName]
			altTmp = dna.BaseToString(dna.ToUpper(currentSeq[snps[i].RefPos-1])) + snps[i].QuerySub
			infoTag = "POS=" + strconv.FormatInt(snps[i].QueryPos, 10)
			curr = &Vcf{Chr: snps[i].RefName, Pos: snps[i].RefPos, Id: snps[i].QueryName, Ref: dna.BaseToString(dna.ToUpper(currentSeq[snps[i].RefPos-1])), Alt: altTmp, Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=INS"}
			answer = append(answer, curr)
		//logic for deletions in vcf record
		case strings.Compare(snps[i].QuerySub, "-") == 0:
			var refTmp string
			currentSeq = reference[snps[i].RefName]
			contigPointer = contigs[snps[i].QueryName]
			refTmp = snps[i].RefSub + dna.BaseToString(dna.ToUpper(contigPointer[snps[i].QueryPos-1]))
			infoTag = "POS=" + strconv.FormatInt(snps[i].QueryPos, 10)
			curr = &Vcf{Chr: snps[i].RefName, Pos: snps[i].RefPos, Id: snps[i].QueryName, Ref: refTmp, Alt: dna.BaseToString(dna.ToUpper(currentSeq[snps[i].RefPos-1])), Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=DEL"}
			answer = append(answer, curr)
		}
	}
	return answer
}

func SnpTruth(s *Snp) int {
	if strings.Compare(s.RefSub, s.QuerySub) != 0 && strings.Compare(s.RefSub, "-") != 0 && strings.Compare(s.QuerySub, "-") != 0 {
		return 0
	} else {
		return -1
	}
}

func MummerToSNP(snpsFile string) ([]*Snp, error) {
	var answer []*Snp
	var curr *Snp
	var line string
	file, err := os.Open(snpsFile)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil, err
	}
	var err2 error
	var rline []byte
	for ; err2 != io.EOF; rline, _, err2 = reader.ReadLine() {
		line = string(rline[:])
		data := strings.Split(line, "\t")
		//fmt.Println("there is data here")
		switch {
		case len(data) == 1:
			//these lines are sequences, and we are not recording them
			//fmt.Println("found sequences")
		case len(line) == 0:
			//blank line
			//fmt.Println("found blank")
		case len(data) == 12:
			//fmt.Println("found header line")

			pRef, _ := strconv.ParseInt(data[0], 10, 64)
			pQuery, _ := strconv.ParseInt(data[3], 10, 64)
			
			if len(answer) == 0 {
				curr = &Snp{RefName: data[10], RefPos: pRef, RefSub: data[1], QueryName: data[11], QueryPos: pQuery, QuerySub: data[2]}
			}

			curr = &Snp{RefName: data[10], RefPos: pRef, RefSub: data[1], QueryName: data[11], QueryPos: pQuery, QuerySub: data[2]}
			answer = append(answer, curr)

		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer, nil

}

/*
func CompareCoord(alpha *Snp, beta *Snp) int {
	if alpha.RefPos < beta.RefPos {
		return -1
	}
	if alpha.RefPos > beta.RefPos {
		return 1
	}
	return 0
}

func CompareName(alpha string, beta string) int {
	return strings.Compare(alpha, beta)
}

func CompareSnp(alpha *Snp, beta *Snp) int {
	compareStorage := CompareName(alpha.RefName, beta.RefName)
	if compareStorage != 0 {
		return compareStorage
	} else {
		return CompareCoord(alpha, beta)
	}
}

func Sort(snps []*Snp) {
	sort.Slice(snps, func(i, j int) bool { return CompareSnp(snps[i], snps[j]) == -1 })
}
*/

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

func Sort(vcfFile []*Vcf) {
	sort.Slice(vcfFile, func(i, j int) bool { return CompareVcf(vcfFile[i], vcfFile[j]) == -1 })
}
//visualization scripts
func PrintVcf(input []*Vcf) {
	for i := range input {
		fmt.Printf("%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", input[i].Chr, input[i].Pos, input[i].Id, input[i].Ref, input[i].Alt, input[i].Qual, input[i].Filter, input[i].Info, input[i].Format, input[i].Unknown)
	}
}

func PrintSnp(input []*Snp) {
	for i := range input {
		fmt.Printf("%s\t%s\t%v\t%v\t%s\t%s\n", input[i].RefName, input[i].QueryName, input[i].RefPos, input[i].QueryPos, input[i].RefSub, input[i].QuerySub)
	}
}



func PrintHeader() {
	var header []string
	t := time.Now()
	header = append(header, "##fileformat=VCFv4.2\n" +
		"##fileDate=" + t.Format("20060102") + "\n" +
		"##source=github.com/vertgenlab/gonomics\n" +
		"##reference=gasAcu1")



	header = append(header, "##contig=<ID=chrI,length=28185914>")
	header = append(header, "##contig=<ID=chrII,length=23295652>")
	header = append(header, "##contig=<ID=chrIII,length=16798506>")
	header = append(header, "##contig=<ID=chrIV,length=32632948>")
	header = append(header, "##contig=<ID=chrIX,length=20249479>")
	header = append(header, "##contig=<ID=chrUn,length=62550211>")
	header = append(header, "##contig=<ID=chrV,length=12251397>")
	header = append(header, "##contig=<ID=chrVI,length=17083675>")
	header = append(header, "##contig=<ID=chrVII,length=27937443>")
	header = append(header, "##contig=<ID=chrVIII,length=19368704>")
	header = append(header, "##contig=<ID=chrX,length=15657440>")
	header = append(header, "##contig=<ID=chrXI,length=16706052>")
	header = append(header, "##contig=<ID=chrXII,length=18401067>")
	header = append(header, "##contig=<ID=chrXIII,length=20083130>")
	header = append(header, "##contig=<ID=chrXIV,length=15246461>")
	header = append(header, "##contig=<ID=chrXIX,length=20240660>")
	header = append(header, "##contig=<ID=chrXV,length=16198764>")
	header = append(header, "##contig=<ID=chrXVI,length=18115788>")
	header = append(header, "##contig=<ID=chrXVII,length=14603141>")
	header = append(header, "##contig=<ID=chrXVIII,length=16282716>")
	header = append(header, "##contig=<ID=chrXX,length=19732071>")
	header = append(header, "##contig=<ID=chrXXI,length=11717487>")
	header = append(header, "##contig=<ID=chrM,length=15742>")
	header = append(header, "##phasing=none\n" +
		"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">\n" + 
		//"##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">\n" +
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
	for i := range header {
		fmt.Println(header[i])
	}

}
