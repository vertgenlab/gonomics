package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"strings"
)

type GgToolsExe struct {
	Cmd     *flag.FlagSet
	Axtfile string
	Vcfs    string
	Out     string
}

func GgtoolsArgs() *GgToolsExe {
	ggT := &GgToolsExe{Cmd: flag.NewFlagSet("ggtools", flag.ExitOnError)}
	ggT.Cmd.StringVar(&ggT.Axtfile, "axt", "", "axt pairwise alignment file used to create Vcfs")
	ggT.Cmd.StringVar(&ggT.Vcfs, "vcf", "", "vcf file combined with fasta reference to make a genome graph")
	ggT.Cmd.StringVar(&ggT.Out, "out", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")

	ggT.Cmd.StringVar(&ggT.Axtfile, "a", "", "axt pairwise alignment file used to create Vcfs")
	ggT.Cmd.StringVar(&ggT.Vcfs, "v", "", "vcf file combined with fasta reference to make a genome graph")
	ggT.Cmd.StringVar(&ggT.Out, "o", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")

	ggT.Cmd.Usage = ggtoolsExtend
	return ggT
}
func RunGgTools() error {
	ggT := GgtoolsArgs()
	Init(ggT.Cmd, os.Args[2:])
	tail := ggT.Cmd.Args()
	graphTools(ggT.Out, ggT.Axtfile, ggT.Vcfs, tail)
	return nil
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
