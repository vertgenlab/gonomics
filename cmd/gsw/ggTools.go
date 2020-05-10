package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"strings"
)

type GgToolsSettings struct {
	Cmd     *flag.FlagSet
	Axtfile string
	Vcfs    string
	Out     string
}

func ggToolsUsage() {
	fmt.Printf(
		"  ggtools\tGenomic utilities to create, manipulate and operate on genome graphs\n")
}

func ggtoolsExtend() {
	fmt.Print(
		"Usage:\n" +
			"  gsw ggtools [options] ref\n\n" +
			"Options:\n" +
			"  -v, --vcf\t\tProvide a VCF to create a graph reference (.gg) used in gsw align\n" +
			"  -a, --axt\t\tUse axt generated from UCSC kentUtils to create a VCF\n\n")
}

func initGgtoolsArgs() *GgToolsSettings {
	ggT := &GgToolsSettings{Cmd: flag.NewFlagSet("ggtools", flag.ExitOnError)}
	ggT.Cmd.StringVar(&ggT.Axtfile, "axt", "", "axt pairwise alignment file used to create Vcfs")
	ggT.Cmd.StringVar(&ggT.Vcfs, "vcf", "", "vcf file combined with fasta reference to make a genome graph")
	ggT.Cmd.StringVar(&ggT.Out, "out", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")

	ggT.Cmd.StringVar(&ggT.Axtfile, "a", "", "axt pairwise alignment file used to create Vcfs")
	ggT.Cmd.StringVar(&ggT.Vcfs, "v", "", "vcf file combined with fasta reference to make a genome graph")
	ggT.Cmd.StringVar(&ggT.Out, "o", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")

	ggT.Cmd.Usage = ggtoolsExtend
	return ggT
}
func RunGgTools() {
	ggT := initGgtoolsArgs()
	ggT.Cmd.Parse(os.Args[2:])
	graphTools(ggT.Out, ggT.Axtfile, ggT.Vcfs, ggT.Cmd.Args())
}

func graphTools(out string, axtfile string, vcfCalls string, files []string) {
	fa := fasta.Read(files[0])
	for i := 0; i < len(files); i++ {
		if strings.HasSuffix(axtfile, ".axt") {
			axtfile := axt.Read(axtfile)
			axt.AxtVcfToFile(out, axtfile, fa)
		} else if strings.HasSuffix(vcfCalls, ".vcf") || strings.HasSuffix(files[i], "vcf.gz") {
			vcfs := vcf.Read(vcfCalls)
			gg := simpleGraph.VariantGraph(fa, vcfs)
			simpleGraph.Write(out, gg)
		}
	}
}
