package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func usage() {
	fmt.Print(
		"callVariants - Calls variants from input sam or giraf based on a linear or graph reference.\n" +
			"Usage:\n" +
			" callVariants [options] \n" +
			"\t-i experimental.sam \n" +
			"\t-n normal.sam \n" +
			"\t-lr reference.fasta \n" +
			"\t-o output.vcf \n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func callVariants(linearRef string, graphRef string, expSamples string, normSamples string, outFile string, afThreshold float64, sigThreshold float64, minMapQ int64, memBufferSize int) {
	var ref interface{}
	output := fileio.MustCreate(outFile)
	defer output.Close()

	if linearRef != "" {
		ref = fasta.Read(linearRef)
	} else if graphRef != "" {
		ref = simpleGraph.Read(graphRef)
	}

	answer := alleles.FindNewVariation(ref, expSamples, normSamples, afThreshold, sigThreshold, minMapQ, memBufferSize)

	vcf.WriteHeader(output)
	for vcfRecord := range answer {
		vcf.PrintSingleLine(vcfRecord)
		vcf.WriteVcf(output, vcfRecord)
	}
}

func main() {
	var outFile *string = flag.String("out", "", "Write output to a file [.vcf].")
	var sigThreshold *float64 = flag.Float64("p", 0.05, "Do not output variants with p value greater than this value.")
	var afThreshold *float64 = flag.Float64("af", 0.01, "Variants with allele frequency less than this value will be treated as p = 1.")
	var linearReference *string = flag.String("lr", "", "Linear reference used for alignment [.fasta].")
	var graphReference *string = flag.String("gr", "", "Graph reference used for alignment [.gg].")
	var experimentalSamples *string = flag.String("i", "", "Input experimental sample(s) [.sam, .giraf]. Can be a file or directory.")
	var normalSamples *string = flag.String("n", "", "Input normal sample(s) [.sam, .giraf]. Can be a file or directory. If no normal samples are given, each experimental sample will me measured against the other experimental samples.")
	var minMapQ *int64 = flag.Int64("minMapQ", 20, "Exclude all reads with mapping quality less than this value")
	var memBufferSize *int = flag.Int("memBuffer", 100, "Maximum number of allele records to store in memory at once")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *linearReference != "" && *graphReference != "" {
		log.Fatalln("ERROR: Cannot input both linear and graph references")
	}

	if *experimentalSamples == "" || *outFile == "" || (*linearReference == "" && *graphReference == "") {
		flag.Usage()
		log.Fatalf("ERROR: Must include parameters for -i, -out, (-lr or -gr)")
	}
	flag.Parse()

	callVariants(*linearReference, *graphReference, *experimentalSamples, *normalSamples, *outFile, *afThreshold, *sigThreshold, *minMapQ, *memBufferSize)
}
