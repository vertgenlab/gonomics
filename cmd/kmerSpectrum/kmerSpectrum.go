package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"sort"
)

type Settings struct {
	FaFile     string
	BedFile    string
	OutFile    string
	KmerLength int
}

func kmerSpectrum(s Settings) {
	var err error
	var currStart int
	var foundInMap bool
	var currKmerBases []dna.Base = make([]dna.Base, s.KmerLength)
	var currKmer, currRevCompKmer, currKey string
	var kmerMap map[string]int
	var currKeys []string
	records := fasta.Read(s.FaFile)
	faMap := fasta.ToMap(records)
	regions := bed.Read(s.BedFile)
	out := fileio.EasyCreate(s.OutFile)
	_, err = fmt.Fprintf(out, "Filename\tRegion\tKmer\tCount\n")
	exception.PanicOnErr(err)
	for currRegion := range regions {
		kmerMap = make(map[string]int, 0)
		if _, foundInMap = faMap[regions[currRegion].Chrom]; !foundInMap {
			log.Fatalf("Error: Chrom in bed region: %v, not found in reference genome.\n", regions[currRegion].Chrom)
		}

		for currStart = regions[currRegion].ChromStart; currStart < regions[currRegion].ChromEnd-s.KmerLength; currStart++ {
			copy(currKmerBases, faMap[regions[currRegion].Chrom][currStart:currStart+s.KmerLength])
			currKmer = dna.BasesToString(currKmerBases)
			dna.ReverseComplement(currKmerBases)
			currRevCompKmer = dna.BasesToString(currKmerBases)
			if _, foundInMap = kmerMap[currKmer]; !foundInMap {
				kmerMap[currKmer] = 1
			} else {
				kmerMap[currKmer]++
			}
			if _, foundInMap = kmerMap[currRevCompKmer]; !foundInMap {
				kmerMap[currRevCompKmer]++
			} else {
				kmerMap[currRevCompKmer]++
			}
		}
		// lexicographic key sort kmerMap to ensure deterministic outfile order
		currKeys = make([]string, 0)
		for currKey = range kmerMap {
			currKeys = append(currKeys, currKey)
		}
		sort.Strings(currKeys)
		for _, currKey = range currKeys {
			_, err = fmt.Fprintf(out, "%s\t%v\t%s\t%v\n", s.BedFile, regions[currRegion].Name, currKey, kmerMap[currKey])
			exception.PanicOnErr(err)
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"kmerSpectrum - Generate k-mer spectra for genomic regions.\n" +
			"Usage:\n" +
			"kmerSpectrum input.fa regions.bed out.txt" +
			"Output txt file will be a tab-delimited table displaying the input filename, region names, kmers, and kmer counts.\n" +
			"As such, the input bed files are expected to have uniquely specified names in the fourth column.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var kmerLength *int = flag.Int("kmerLength", 3, "Set the length of k-mers. Must be non-negative.\n")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	faFile := flag.Arg(0)
	bedFile := flag.Arg(1)
	outFile := flag.Arg(2)

	s := Settings{
		FaFile:     faFile,
		BedFile:    bedFile,
		OutFile:    outFile,
		KmerLength: *kmerLength,
	}
	kmerSpectrum(s)
}
