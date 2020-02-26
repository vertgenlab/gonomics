package main

//TODO: still a work in progress for the gsw aligning executable
import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"

	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"gsw - A genomics toolkit that operates on sequence graphs.\n" +
			"\nusage:\n" +
			"\tgsw [options] ref.fa/.gg\n" +
			"options:\n" +
			"\t--align input.fastq.gz\n" +
			"\t\tgraph-smith-waterman alignment algorithm: align sequencing\n\t\tdata to a sequence graph reference\n" +
			"\t\tmust have ref.gg or create a genome graph\n" +
			"\t\tusing ref.fa and .axt or .vcf\n" +
			"\t--vcf input.vcf\n\t\tuse vcf as input to split reference\n\t\tinto a sequence graph with nodes.\n" +
			"\t--axt input.axt\n" +
			"\t\tuse axt pairwise alignment to convert to vcf,\n\t\tthen split reference into a sequence graph with nodes\n" +
			"\t--out output[.sam/.vcf/.gg]\n" +
			"\t\t.sam format output from graph alignment\n" +
			"\t\t.vcf small snp/ins/del from axt to vcf\n" +
			"\t\t.gg write genome graph to file\n" +
			"\t--test number\n" +
			"\t\tSimulate n number of reads from given fasta and realigns reads\n")
	//flag.PrintDefaults()
}

//WORK IN PROGRESS
func main() {
	var vcfs []*vcf.Vcf
	var fa []*fasta.Fasta
	var gg *simpleGraph.SimpleGraph
	var tagAxt *string = flag.String("axt", "", "axt alignment file")
	var vcfTag *string = flag.String("vcf", "", "vcf file")
	var outTag *string = flag.String("out", "", "final output, .vcf/.gg/.sam")
	var alignFlag *bool = flag.Bool("align", false, "in.fastq out.sam")
	//var fqFlag *string = flag.String("-fastq", "", "read.fastq")
	var testSim *int = flag.Int("test", 0, "Simulates n reads from and given fasta reference")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	ref := flag.Arg(0)
	if ref == "" {
		usage()
	}
	if *tagAxt != "" {
		vcfs = mkVcf(*tagAxt)
		if *outTag != "" {
			if strings.HasSuffix(*outTag, ".vcf") {
				vcf.Write(*outTag, vcfs)
			}
		}
		if *outTag == "" {
			vcf.PrintVcf(vcfs)
		}
	}
	if *vcfTag != "" {
		vcfs = vcf.Read(*vcfTag)
		if strings.HasSuffix(ref, ".fa") || strings.HasPrefix(ref, ".fasta") {
			fa = fasta.Read(ref)
			gg = simpleGraph.FaToGenomeGraph(fa, vcfs)
			if strings.HasSuffix(*outTag, ".gg") {
				simpleGraph.Write(*outTag, gg)
			}
			if *outTag == "" {
				simpleGraph.PrintGraph(gg)
			}
		}
	}
	if *testSim > 0 {
		gg = simpleGraph.Read("gasAcu1.fa")

		simpleGraph.WriteReadsToFile(gg.Nodes, 150, 10000)
		//simpleGraph.GenomeDiversitySimulator()

	}
	if *alignFlag == true {
		if strings.HasSuffix(ref, ".gg") {
			gg = simpleGraph.Read(ref)
		}
		if strings.HasSuffix(ref, ".fa") {
			gg = simpleGraph.Read(ref)
		}
		//flag args(0) is usually the reference
		reads := flag.Args()
		//reads = reads[1:]
		//if user only provides single fastq, aligns single read
		if len(reads) == 2 {
			//if usuer provides a .sam as output will write alignment
			//to that file; otherwise, software will write to STDout
			if strings.HasSuffix(*outTag, ".sam") {
				simpleGraph.GSWsBatchDraft(gg, reads[1], *outTag)
			}
		}
		if len(reads) == 3 {
			if strings.HasSuffix(*outTag, ".sam") {
				simpleGraph.GSWsBatchPair(gg, reads[1], reads[2], *outTag)
			}
		}
		//if user provides paired end reads, will do alignment on paired end reads
	}
}

func mkVcf(filename string) []*vcf.Vcf {
	axtFile := axt.Read(filename)
	vcfs := axt.CallSnpsToVcf(axtFile)
	vcf.Sort(vcfs)
	return vcfs
}
