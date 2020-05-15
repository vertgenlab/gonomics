package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/vcf"
)

type globalInit struct {
	Cmd        *flag.FlagSet
	NonOverlap bool
	Extend 	   int
	Format     string
	Threads    int
	Output     string
	Help       string
}

type regions struct {
	Chr    string
	Start  int
	End    int
	Record *Data
}

type Data struct {
	axts *axt.Axt
	beds *bed.Bed
	vcfs *vcf.Vcf
}

func settings(name string) *globalInit {
	overlapGenome := &globalInit{Cmd: flag.NewFlagSet(name, flag.ExitOnError)}
	overlapGenome.Cmd.Usage = usage

	overlapGenome.Cmd.BoolVar(&overlapGenome.NonOverlap, "nonOverlap", false, "return non-overlapping regions")
	overlapGenome.Cmd.BoolVar(&overlapGenome.NonOverlap, "n", false, "return non-overlapping regions")

	overlapGenome.Cmd.IntVar(&overlapGenome.Extend, "extend", 0, "add to start and end of target coordinates, requires a chrom.sizes file")
	overlapGenome.Cmd.IntVar(&overlapGenome.Extend, "e", 0, "add to start and end of target coordinates, requires a chrom.sizes file")

	overlapGenome.Cmd.StringVar(&overlapGenome.Format, "format", "", "choose between your Target fmt or Bed as final output")
	overlapGenome.Cmd.StringVar(&overlapGenome.Format, "f", "", "choose between your Target fmt or Bed as final output")

	overlapGenome.Cmd.IntVar(&overlapGenome.Threads, "threads", 1, "number of CPUs for Goroutine concurrency")
	overlapGenome.Cmd.IntVar(&overlapGenome.Threads, "t", 1, "number of CPUs for Goroutine concurrency")

	overlapGenome.Cmd.StringVar(&overlapGenome.Output, "out", "/dev/stdout", "filename of final data")
	overlapGenome.Cmd.StringVar(&overlapGenome.Output, "o", "/dev/stdout", "filename of final data")

	overlapGenome.Cmd.StringVar(&overlapGenome.Help, "help", "", "view detailed help message specified option")
	overlapGenome.Cmd.StringVar(&overlapGenome.Help, "h", "", "view detailed help message specified option")
	return overlapGenome
}
