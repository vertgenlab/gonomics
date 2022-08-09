package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func haplotypeGenerator(genomeFile string, snpFile string, regionFile string, outdir string) {
	variants, header := vcf.Read(snpFile)
	variantIntervals := make([]interval.Interval, len(variants))
	for i := range variants {
		variantIntervals[i] = variants[i]
	}

	tree := interval.BuildTree(variantIntervals)

	regions := bed.Read(regionFile)
	genome := fasta.Read(genomeFile)

	for i := range regions {
		overlappingVariants := interval.Query(tree, regions[i], "any")
		refHaplotype := fasta.Extract(genome, regions[i].ChromStart, regions[i].ChromEnd, regions[i].Chrom)
		outputFilename := fmt.Sprintf("%s/%s.%v.%v.fa", outdir, regions[i].Chrom, regions[i].ChromStart, regions[i].ChromEnd)
		out := fileio.EasyCreate(outputFilename)
		sampleHaplotypes := make([]fasta.Fasta, len(header.Samples)*2)
		for j := range sampleHaplotypes {
			copy(sampleHaplotypes[j], refHaplotype)
			for k := range overlappingVariants {
				// For each variant, update corresponding position in each sample haplotype (position = position - ChromStart - 1)

			}
			fasta.WriteFasta(out, sampleHaplotypes[j], 50)
		}

		err := out.Close()
		exception.PanicOnErr(err)

	}

}

func usage() {
	fmt.Print(
		"haplotypeGenerator - Generate unique haplotypes for provided regions from genetic variation data.\n" +
			"Ignores structural variants, considers only substitutions\n" +
			"Usage:\n" +
			"genomeFile.fa snpFile.vcf regionFile.bed outDir\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	genomeFile := flag.Arg(0)
	snpFile := flag.Arg(1)
	regionFile := flag.Arg(2)
	outdir := flag.Arg(3)

	haplotypeGenerator(genomeFile, snpFile, regionFile, outdir)
}
