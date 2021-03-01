package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"strings"
)

type GgToolsSettings struct {
	Cmd       *flag.FlagSet
	TargetFa  string
	QueryFa   string
	FmtOutput string
	Out       string
}

func ggToolsUsage() {
	fmt.Printf(
		"  ggtools\tGenomic utilities to create, manipulate and operate on genome graphs\n")
}

func ggtoolsExtend() {
	fmt.Print(
		"Usage:\n" +
			"  gsw ggtools [options] input.file\n\n" +
			"Required:\n" +
			"  -t, --target\tTarget reference fasta file\n" +
			"  -q, --query\tQuery fasta file required for chain file inputs\n\n" +
			"Options:\n" +
			"  -f, --format\t  Pick file format for output, only applies to file conversions utils\n" +
			"  -o, --out\t  [.chain/.gg/.sam/.vcf/.gz]  (default: /dev/stdout)\n\n")
}

func initGgtoolsArgs() *GgToolsSettings {
	ggT := &GgToolsSettings{Cmd: flag.NewFlagSet("ggtools", flag.ExitOnError)}
	ggT.Cmd.StringVar(&ggT.FmtOutput, "format", "", "Pick file format for output, only applies to file conversions utils")
	ggT.Cmd.StringVar(&ggT.TargetFa, "target", "", "Specify target reference fasta file")
	ggT.Cmd.StringVar(&ggT.QueryFa, "query", "", "Specify query fasta file")
	ggT.Cmd.StringVar(&ggT.Out, "out", "/dev/stdout", "Output filename, [.gg/.vcf]")

	ggT.Cmd.StringVar(&ggT.TargetFa, "t", "", "Specify target reference fasta file")
	ggT.Cmd.StringVar(&ggT.QueryFa, "q", "", "Specify query fasta file")
	ggT.Cmd.StringVar(&ggT.FmtOutput, "f", "", "Pick file format for output, only applies to file conversions utils")
	ggT.Cmd.StringVar(&ggT.Out, "o", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")
	ggT.Cmd.Usage = ggtoolsExtend
	return ggT
}
func RunGgTools() {
	ggT := initGgtoolsArgs()
	ggT.Cmd.Parse(os.Args[2:])
	if len(ggT.Cmd.Args()) != 1 {
		ggT.Cmd.Usage()
	} else {
		inFile := ggT.Cmd.Arg(0)
		log.Printf("Reading: %s\n", inFile)
		switch true {
		case chain.IsChainFile(inFile):
			log.Printf("Input chain detected...\n")
			if strings.Compare(ggT.TargetFa, "") != 0 && strings.Compare(ggT.QueryFa, "") != 0 {
				log.Printf("Converting chain to %s...\n", ggT.FmtOutput)
				convertChains(inFile, ggT.TargetFa, ggT.QueryFa, ggT.FmtOutput, ggT.Out)
			} else {
				ggT.Cmd.Usage()
				log.Fatalf("Error: Must specify both target and query fasta files...\n")
			}
		case vcf.IsVcfFile(inFile):
			log.Printf("Input vcf detected...\n")
			if strings.Compare(ggT.TargetFa, "") != 0 {
				log.Printf("Converting vcf to %s...\n", ggT.FmtOutput)
				simpleGraph.Write(ggT.Out, vcfToSimpleGraph(inFile, ggT.TargetFa))
			} else {
				ggT.Cmd.Usage()
				log.Fatalf("Error: Must specify target reference fasta file...\n")
			}
		case axt.IsAxtFile(inFile):
			log.Printf("Input axt detected...\n")
			if strings.Compare(ggT.TargetFa, "") != 0 {
				log.Printf("Converting axt alignment to %s...\n", ggT.FmtOutput)
				convertAxt(inFile, ggT.FmtOutput, ggT.TargetFa, ggT.Out)
			}
		default:
			// /ggT.Cmd.Usage()
			errorMessage()
		}
		log.Printf("Finished writing %s file\n-xoxo gg\n", ggT.FmtOutput)
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

func chainToBed(chainFmt, bedFmt string, target bool, out string) {
	reader, _ := chain.GoReadToChan(chainFmt)

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
		if !strings.Contains(data.Ref, "N") && !strings.Contains(data.Alt[0], "N") {
			chrMap[data.Chr] = append(chrMap[data.Chr], data)
		}
	}
	return chrMap
}

func FaVcfChannels(ref string, vcfInput string) *simpleGraph.SimpleGraph {
	vcfFilteredMap := vcfSplitChrNoN(vcfInput)
	faReader := fasta.GoReadToChan(ref)
	gg := simpleGraph.VariantGraph(faReader, vcfFilteredMap)
	return gg
}

/*
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
/*
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
}*/
