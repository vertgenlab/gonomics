package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Print(
		"callVariants - Calls variants from input sam or giraf based on a linear or graph reference.\n" +
			"Usage:\n" +
			" callVariants [options] \n" +
			"\t-i experimental.sam \n" +
			"\t-n normal.sam \n" +
			"\t-r reference.fasta \n" +
			"\t-o output.vcf \n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var outFile *string = flag.String("out", "stdout", "Write output to a file [.vcf].")
	var sigThreshold *float64 = flag.Float64("p", 0.05, "Do not output variants with p value greater than this value.")
	var afThreshold *float64 = flag.Float64("af", 0.01, "Do not output variants with allele frequency less than this value.")
	var reference *string = flag.String("r", "", "Reference used for alignment [.fasta, .gg]. Can be linear or graph.")
	var experimentalSamples *string = flag.String("i", "", "Input experimental sample(s) [.sam, .giraf]. Can be a file or directory.")
	var normalSamples *string = flag.String("n", "", "Input normal sample(s) [.sam, .giraf]. Can be a file or directory. If no normal samples are given, each experimental sample will me measured against the other experimental samples.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *experimentalSamples == "" || *reference == "" || *outFile == "" {
		flag.Usage()
		log.Fatalf("ERROR: Must include parameters for -i, -r, -o")
	}
	flag.Parse()

	fmt.Sprintf("test", normalSamples, experimentalSamples, reference, afThreshold, sigThreshold, outFile)
}
