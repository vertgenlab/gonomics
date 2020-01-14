package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	//"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/fastq"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"gsw - A genomics toolkit that operates on sequence graphs.\n" +
			"\nusage:\n" +
			"\tgsw ref.fa/.gg [options] --vcf/axt file\n" +
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
			"\t\t.gg write genome graph to file\n")
	flag.PrintDefaults()
}

func main() {

	//var i int
	//var tileSize int = 12
	//var stepSize int = 1

	//var mappedRead *sam.SamAln
	var vcfs []*vcf.Vcf
	var fa []*fasta.Fasta

	var gg *simpleGraph.SimpleGraph

	var tagAxt *string = flag.String("axt", "", "axt alignment file")
	var vcfTag *string = flag.String("vcf", "", "vcf file")
	var outTag *string = flag.String("out", "", "final output, .vcf/.gg/.sam")
	var alignFlag *string = flag.String("align", "", "in.fastq out.sam")
	//var fqFlag *string = flag.String("-fastq", "", "read.fastq")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	ref := flag.Arg(0)

	/*
		if *tagAxt == "" && *vcfTag == "" && *outTag == "" && *alignFlag == "" && *fqFlag == "" {
			flag.Usage()
		}*/

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
	if *alignFlag != "" {
		if strings.HasSuffix(ref, ".gg") {
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
				simpleGraph.GSWAlignerOne(gg, reads[1], *outTag)
			}
		}
		//if user provides paired end reads, will do alignment on paired end reads

	}

	/*
		var gg *simpleGraph.SimpleGraph

		if *vcfTag != "" {
			vcfs = vcf.Read(*vcfTag)
			if strings.HasSuffix(ref, ".fa") || strings.HasPrefix(ref, ".fasta") {
				fa := fasta.Read(ref)
				gg = simpleGraph.FaToGenomeGraph(fa, vcfs)
				//simpleGraph.PrintGraph(gg)
				if strings.HasSuffix(ref, ".gg") {
					simpleGraph.Write(*outTag, gg)
				}
			}
		}
		if *alignFlag != "" {
			if *tagAxt != "" {
				fa := fasta.Read(ref)
				gg = simpleGraph.FaToGenomeGraph(fa, vcfs)
			}

			//&& *fqFlag != "" && strings.HasSuffix(ref, ".sam")
			log.Printf("Indexing the genome graph...\n")
			ham5 := simpleGraph.IndexGenomeDev(gg.Nodes, tileSize, stepSize)
			if *fqFlag == "" {
				log.Fatal("alignment setting set, but no fastq file was detected")
			} else {
				fq := fastq.Read(*fqFlag)
				m, trace := simpleGraph.SwMatrixSetup(10000)
				for i = 0; i < len(fq);i++ {
					mappedRead = simpleGraph.GraphSmithWaterman(gg, fq[i], ham5, tileSize, m, trace)
					fmt.Printf("%s\n", sam.SamAlnToString(mappedRead))
				}
			}
			/*
			log.Printf("Aligning reads...\n")
			var done bool
			var fq *fastq.Fastq

			outFile, _ := os.Create(*outTag)
			defer outFile.Close()
			sam.WriteHeaderToFileHandle(outFile, header)
			file := fileio.EasyOpen(*fqFlag)
			defer file.Close()

			mSet, trace := swMatrixSetup(10000)

			for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {

			}*/

	//}
}

func logicTags(tagAxt string, outTag string) {

}

func mkVcf(filename string) []*vcf.Vcf {
	axtFile := axt.Read(filename)
	vcfs := axt.CallSnpsToVcf(axtFile)
	vcf.Sort(vcfs)
	return vcfs
}
