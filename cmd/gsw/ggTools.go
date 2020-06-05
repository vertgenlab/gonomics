package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	//"log"
	"strings"
	"sync"
)

type GgToolsSettings struct {
	Cmd           *flag.FlagSet
	Axtfile       string
	BedFmt        string
	ChainFmt      string
	VcfFmt        string
	NonOverlap    bool
	SelectOverlap bool
	TargetChrom   string
	QueryChrom    string
	FmtOutput     string
	Out           string
}

func ggToolsUsage() {
	fmt.Printf(
		"  ggtools\tGenomic utilities to create, manipulate and operate on genome graphs\n")
}

func ggtoolsExtend() {
	fmt.Print(
		"Usage:\n" +
			"  gsw ggtools [options] ref\n\n" +
			"Required:\n" +
			"  Fasta reference\n\n" +
			"Options:\n" +
			"  -a, --axt\tConvert axt generated from UCSC kentUtils to a different file format\n\t\tExample:\tgsw ggtools -a pairwise.axt -f vcf -o genome.vcf.gz ref.fa\n\t\tOutfmt:\t\t[.vcf/.gz/.chain/.bed]\n" +
			"  -b  --bed\tAnalyze bed formated genomic regions\n\t\tComing Soon!\n" +
			"  -c, --chain\tInput a chain file for addional graph and assembly capabilities\n\t\tComing Soon!\n" +
			"  -tLen  --targetSizes\tA two column-tab file target name and length, required for converting axt to chain\n" +
			"  -qLen  --querySizes\t" +
			"  -v, --vcf\tProvide a VCF to create a graph reference (.gg) used in gsw align\n\t\tExample:\tgsw ggtools -v variance.vcf.gz -o variance.gg ref.fa\n" +
			"Settings:\n" +
			"  -s, --select\t  Query bed regions that overlap target axt or chain formats\n" +
			"  -i, --invert\t  Filter nonOverlapping bed regions from target axt or chain formats\n" +
			"  -f, --format\t  Pick file format for output, only applies to file conversions utils\n" +
			"  -o, --out\t  Filename[.bed/.chain/.gg/.sam/.vcf/.gz/]  (default: /dev/stdout)\n\n")
}

func initGgtoolsArgs() *GgToolsSettings {
	ggT := &GgToolsSettings{Cmd: flag.NewFlagSet("ggtools", flag.ExitOnError)}
	ggT.Cmd.StringVar(&ggT.Axtfile, "axt", "", "Convert axt generated from UCSC kentUtils to a different file format")
	ggT.Cmd.StringVar(&ggT.BedFmt, "bed", "", "Analyze bed formated genomic regions")
	ggT.Cmd.StringVar(&ggT.ChainFmt, "chain", "", "Input a chain file for addional graph and assembly capabilities")
	ggT.Cmd.StringVar(&ggT.FmtOutput, "format", "", "Pick file format for output, only applies to file conversions utils")
	ggT.Cmd.StringVar(&ggT.TargetChrom, "targetSizes", "", "A two column-tab file containing target name and length is required for converting axt to chain")
	ggT.Cmd.StringVar(&ggT.QueryChrom, "querySizes", "", "A two column-tab file containing query name and length is required for converting axt to chain")

	ggT.Cmd.BoolVar(&ggT.SelectOverlap, "select", false, "Query bed regions that overlap target axt or chain formats")
	ggT.Cmd.BoolVar(&ggT.NonOverlap, "invert", false, "Filter nonOverlapping bed regions from target axt or chain formats")
	ggT.Cmd.StringVar(&ggT.VcfFmt, "vcf", "", "vcf file combined with fasta reference to make a genome graph")
	ggT.Cmd.StringVar(&ggT.Out, "out", "/dev/stdout", "Output filename, [.gg/.vcf]")

	ggT.Cmd.StringVar(&ggT.Axtfile, "a", "", "Convert axt generated from UCSC kentUtils to a different file format")
	ggT.Cmd.StringVar(&ggT.BedFmt, "b", "", "Analyze bed formated genomic regions")
	ggT.Cmd.StringVar(&ggT.ChainFmt, "c", "", "Input a chain file for addional graph and assembly capabilities")
	ggT.Cmd.StringVar(&ggT.FmtOutput, "f", "", "Pick file format for output, only applies to file conversions utils")
	ggT.Cmd.StringVar(&ggT.TargetChrom, "tLen", "", "A two column-tab file containing target name and length is required for converting axt to chain")
	ggT.Cmd.StringVar(&ggT.QueryChrom, "qLen", "", "A two column-tab file containing query name and length is required for converting axt to chain")

	ggT.Cmd.BoolVar(&ggT.SelectOverlap, "s", false, "Query bed regions that overlap target axt or chain formats")
	ggT.Cmd.BoolVar(&ggT.NonOverlap, "i", false, "Query nonoverlapping bed regions from target axt or chain formats")
	ggT.Cmd.StringVar(&ggT.VcfFmt, "v", "", "vcf file combined with fasta reference to make a genome graph")
	ggT.Cmd.StringVar(&ggT.Out, "o", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")

	ggT.Cmd.Usage = ggtoolsExtend
	return ggT
}
func RunGgTools() {
	ggT := initGgtoolsArgs()
	ggT.Cmd.Parse(os.Args[2:])
	if len(os.Args) < 2 {
		ggT.Cmd.Usage()
	} else {
		graphTools(ggT.Out, ggT.Axtfile, ggT.VcfFmt, ggT.FmtOutput, ggT.Cmd.Arg(0), ggT.BedFmt, ggT.ChainFmt, ggT.SelectOverlap, ggT.NonOverlap, ggT.TargetChrom, ggT.QueryChrom)
	}
}

//This is similar to the function we usually run in main() which takes the flag pointers as arguments and runs our analysis
func graphTools(out string, axtfile string, vcfCalls string, outfmt string, faRef string, bedFmt string, chainFmt string, overlap bool, invert bool, targetChrom string, queryChrom string) {
	switch true {
	case strings.HasSuffix(axtfile, ".axt") && strings.Contains(outfmt, "chain") && strings.Compare(targetChrom, "") != 0 && strings.Compare(queryChrom, "") != 0:
		axtToChain(axtfile, out, targetChrom, queryChrom)
	case overlap && strings.HasSuffix(axtfile, ".axt") && strings.HasSuffix(bedFmt, ".bed"):
		findAxtBedOverlap(axtfile, bedFmt, out, invert)
	case strings.HasSuffix(chainFmt, ".chain") && strings.Contains(outfmt, "bed"):
		//chain.OverlapChainBed(chainFmt)
		//TODO: add feature to select target or query
		chainToBed(chainFmt, outfmt, true, out)
	default:
		errorMessage()
	}
	//These are the functions that require a fasta reference, need to perform check separate from switch case because the value could be empty
	if isFasta(faRef) {
		if strings.HasSuffix(vcfCalls, ".vcf") || strings.HasSuffix(vcfCalls, "vcf.gz") {
			//user can specify as their out format, but if the function made it this far the only option is to create genome graph
			if strings.Contains(outfmt, "graph") || outfmt == "" {
				gg := FaVcfChannels(faRef, vcfCalls)
				simpleGraph.Write(out, gg)
			}
		} else if strings.HasSuffix(axtfile, ".axt") && isFasta(faRef) {
			if strings.Contains(outfmt, "vcf") || outfmt == "" {
				fa := fasta.Read(faRef)
				axtfile := axt.Read(axtfile)
				axt.AxtVcfToFile(out, axtfile, fa)
			}
		} else {
			errorMessage()
		}
	}
}

func isFasta(filename string) bool {
	if strings.HasSuffix(filename, ".fasta") || strings.HasSuffix(filename, ".fa") || strings.HasSuffix(filename, ".fasta.gz") || strings.HasSuffix(filename, ".fa.gz") {
		return true
	} else {
		//log.Fatalf("Error a fasta reference must be provided this is specific file conversion\n")
		return false
	}
}

func axtToChain(axtFmt string, chainFmt string, target string, query string) {
	input := fileio.EasyOpen(axtFmt)

	defer input.Close()

	reader, writer := make(chan *axt.Axt), make(chan *chain.Chain)

	go axt.ReadToChan(input, reader)

	var wg sync.WaitGroup
	wg.Add(1)
	go chain.WriteToFile(chainFmt, writer, nil, &wg)
	var i int = 0
	targetChrom := chromInfo.ReadToMap(target)
	queryChrom := chromInfo.ReadToMap(query)
	for each := range reader {
		writer <- chain.AxtToChain(each, int(targetChrom[each.RName].Size), int(queryChrom[each.QName].Size), i)
		i++
	}
	close(writer)
	wg.Wait()
}

/*
func routineWorker(reader <-chan *chain.Chain, writing chan <- *axt.Axt, target string, query string) {

}*/

func chainToBed(chainFmt, bedFmt string, target bool, out string) {
	input := fileio.EasyOpen(chainFmt)
	defer input.Close()
	reader := make(chan *chain.Chain)
	go chain.ReadToChan(input, reader)

	outFile := fileio.MustCreate(out)
	defer outFile.Close()
	for each := range reader {
		bed.WriteBed(outFile, chain.ChainToBed(each, target), 5)
	}
}

/*
func chainToAxt(chainFmt string, axtFmt string) {
	input := fileio.EasyOpen(axtFmt)
	defer input.Close()
	reader, writer := make(chan *chain.Chain), make(chan *axt.Axt)
	go chain.ReadToChan(input, reader)
	var i int = 0
	var wg sync.WaitGroup
	wg.Add(1)
	for each := range reader {
		writer <- chain.ChainToAxt(each, i)
		i++
	}

}*/

func vcfSplitChrNoN(vcfInput string) map[string][]*vcf.Vcf {
	reader := fileio.EasyOpen(vcfInput)
	defer reader.Close()
	vcf.ReadHeader(reader)
	var data *vcf.Vcf
	var done bool
	chrMap := make(map[string][]*vcf.Vcf)
	for data, done = vcf.NextVcf(reader); !done; data, done = vcf.NextVcf(reader) {
		if !strings.Contains(data.Ref, "N") && !strings.Contains(data.Alt, "N") {
			chrMap[data.Chr] = append(chrMap[data.Chr], data)
		}
	}
	return chrMap
}

func FaVcfChannels(ref string, vcfInput string) *simpleGraph.SimpleGraph {
	vcfFilteredMap := vcfSplitChrNoN(vcfInput)
	faReader := make(chan *fasta.Fasta)
	go fasta.ReadToChan(ref, faReader, nil, true)
	gg := simpleGraph.VariantGraph(faReader, vcfFilteredMap)
	return gg
}

//axt is the select/target regions we are looking for overlaps from
func findAxtBedOverlap(axtFmt string, bedFmt string, output string, invert bool) {
	query := make(chan *bed.Bed)
	go bed.ReadToChan(bedFmt, query)

	target := mkAxtMap(axtFmt)

	answer := fileio.MustCreate(output)
	defer answer.Close()
	for eachRegion := range query {
		if overlapAxtCompareBeds(target[eachRegion.Chrom], eachRegion) {
			//TODO: consider changing this so people can input n fields
			bed.WriteBed(answer, eachRegion, 5)
		} else {
			if invert {
				bed.WriteBed(answer, eachRegion, 5)
			}
		}
	}
}

func mkAxtMap(axtFmt string) map[string][]*axt.Axt {
	input := fileio.EasyOpen(axtFmt)
	reader := make(chan *axt.Axt)

	target := make(map[string][]*axt.Axt)

	go axt.ReadToChan(input, reader)
	for each := range reader {
		target[each.RName] = append(target[each.RName], each)
	}
	return target
}

//check if a single bed record has overlap
func overlapAxtCompareBeds(target []*axt.Axt, query *bed.Bed) bool {
	for _, each := range target {
		if axt.OverlapAxtBed(each, query) {
			return true
		}
	}
	return false
}

//check if a single bed record has overlap
func overlapChainCompareBeds(target []*chain.Chain, query *bed.Bed) bool {
	for _, each := range target {
		//TODO: can only select regions from target at the moment, will add query option
		if chain.OverlapChainBed(each, query, true) {
			return true
		}
	}
	return false
}
