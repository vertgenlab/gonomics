// Command Group: "Statistics & Population Genetics"

// Generate unique haplotypes for provided regions from genetic variation data
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

type Settings struct {
	ReferenceGenomeFile string
	VcfFile             string
	RegionBedFile       string
	IncludeReference    bool
	OutDir              string
	Verbose             int
	LineLength          int
}

func haplotypeGenerator(s Settings) {
	var currState int
	regions := bed.Read(s.RegionBedFile)
	genome := fasta.Read(s.ReferenceGenomeFile)
	genomeMap := helperFastaIndex(genome)
	var currHaplotype fasta.Fasta
	for i := range regions {
		overlappingVariants := make([]vcf.Vcf, 0)
		variantChan, header := vcf.GoReadToChan(s.VcfFile)
		samplesNames := vcf.SampleNamesInOrder(header)
		if s.Verbose > 0 {
			fmt.Printf("working on region index %v.\n", i)
		}
		for currVariant := range variantChan {
			if interval.Overlap(currVariant, regions[i]) && vcf.IsSubstitution(currVariant) {
				overlappingVariants = append(overlappingVariants, currVariant)
			}
		}

		currIndex := genomeMap[regions[i].Chrom]
		refHaplotype := fasta.Extract(genome[currIndex], regions[i].ChromStart, regions[i].ChromEnd, regions[i].Chrom)
		outputFilename := fmt.Sprintf("%s/%s.%v.%v.fa", s.OutDir, regions[i].Chrom, regions[i].ChromStart, regions[i].ChromEnd)
		out := fileio.EasyCreate(outputFilename)

		if s.IncludeReference {
			fasta.WriteFasta(out, refHaplotype, s.LineLength)
		}

		if s.Verbose > 0 {
			fmt.Printf("output filename: %s\n", outputFilename)
			fmt.Printf("last line of header: \n%s\n", header.Text[len(header.Text)-1])
			fmt.Printf("length of header.Samples: %v.\n", len(header.Samples))
		}

		for j := range samplesNames {
			// first we calculate the haplotype for A, and write it out
			currHaplotype = fasta.Copy(refHaplotype)
			currHaplotype.Name = fmt.Sprintf("%s_A", samplesNames[j])
			if s.Verbose > 0 {
				fmt.Printf("working on haplotype: \n%s", currHaplotype.Name)
			}
			for k := range overlappingVariants {
				currState = int(overlappingVariants[k].Samples[j].Alleles[0])
				if currState > 0 {
					currHaplotype.Seq[overlappingVariants[k].Pos-regions[i].ChromStart-1] = dna.StringToBase(overlappingVariants[k].Alt[currState-1])
				}
			}
			fasta.WriteFasta(out, currHaplotype, s.LineLength)

			// now we calculate the haplotype for B, and write it out
			currHaplotype = fasta.Copy(refHaplotype)
			currHaplotype.Name = fmt.Sprintf("%s_B", samplesNames[j])
			if s.Verbose > 0 {
				fmt.Printf("working on haplotype: \n%s", currHaplotype.Name)
			}
			for k := range overlappingVariants {
				currState = int(overlappingVariants[k].Samples[j].Alleles[1])
				if currState > 0 {
					currHaplotype.Seq[overlappingVariants[k].Pos-regions[i].ChromStart-1] = dna.StringToBase(overlappingVariants[k].Alt[currState-1])
				}
			}
			fasta.WriteFasta(out, currHaplotype, s.LineLength)
		}

		err := out.Close()
		exception.PanicOnErr(err)
	}
}

// helperFastaIndex is a helper function that returns a map connecting chromosome names to their index in a fasta slice.
func helperFastaIndex(genome []fasta.Fasta) map[string]int {
	var answer = make(map[string]int)
	for i := range genome {
		answer[genome[i].Name] = i
	}
	return answer
}

func usage() {
	fmt.Print(
		"haplotypeGenerator - Generate unique haplotypes for provided regions from genetic variation data.\n" +
			"Ignores structural variants, considers only substitutions\n" +
			"Usage:\n" +
			"haplotypeGenerator genomeFile.fa snpFile.vcf regionFile.bed outDir\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var lineLength *int = flag.Int("lineLength", 50, "Specify the number of bases per line in the output fasta.\n")
	var includeReference *bool = flag.Bool("includeReference", false, "If set, the program will also write out the reference haplotype as the first result. Note that as only substitutions are considered, this renders the output file a valid multiFa.\n")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 to reveal debug prints.\n")

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
	outDir := flag.Arg(3)

	s := Settings{
		ReferenceGenomeFile: genomeFile,
		VcfFile:             snpFile,
		RegionBedFile:       regionFile,
		OutDir:              outDir,
		IncludeReference:    *includeReference,
		LineLength:          *lineLength,
		Verbose:             *verbose,
	}

	haplotypeGenerator(s)
}
