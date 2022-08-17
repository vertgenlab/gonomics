package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/lastZWriter"
	"io/fs"
	"log"
	"path/filepath"
)

//TODO: write main function with options etc.
func MakeArray(lastZ string, pairwise string, speciesListFile string, refListFile string, allDists string, outText string, outCsv string) {
	speciesList := fileio.EasyOpen(speciesListFile)
	refList := fileio.EasyOpen(refListFile)
	csvFile := fileio.EasyCreate(outCsv)
	fileio.EasyCreate(outText)
	var spec, ref string
	var speciesDone, refDone bool
	var parameters []string
	var matrix string
	var dist int
	var allLines []string
	csvOut := csv.NewWriter(csvFile)
	for ref, refDone = fileio.EasyNextRealLine(refList); !refDone; ref, refDone = fileio.EasyNextRealLine(refList) {
		for spec, speciesDone = fileio.EasyNextRealLine(speciesList); !speciesDone; spec, speciesDone = fileio.EasyNextRealLine(speciesList) {
			if spec != ref {
				parameters, matrix, dist = lastZWriter.AlignSetUp(pairwise, spec, ref, allDists)
				allLines = writeFiles(lastZ, pairwise, ref, spec, parameters, matrix, dist, allLines, csvOut)
			}
		}
	}
	fileio.Write(outText, allLines)
}

func writeFiles(lastZ string, pairwise string, reference string, species string, parameters []string, matrix string, dist int, allLines []string, csvOut *csv.Writer) (lines []string) {
	var currLine string
	par := fmt.Sprintf("%s %s %s %s %s %s %s ", parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], parameters[6])
	allVars := []string{
		reference, species, string(dist), matrix, par}
	csvOut.Write(allVars)

	currLine = fastaFinder(lastZ, pairwise, reference, species, par, matrix)
	allLines = append(allLines, currLine)

	return allLines
}

func fastaFinder(lastZ string, pairwise string, reference string, species string, par string, matrix string) (line string) {
	var currLine string
	tPath := pairwise + "/" + reference
	qPath := pairwise + "/" + species

	filepath.WalkDir(tPath, func(f string, d fs.DirEntry, err error) error {
		if err != nil {
			log.Fatal(err)
		}
		matched, err := filepath.Match("*.fa", f)
		if err != nil {
			return err
		} else if matched {
			tFasta := f

			filepath.WalkDir(qPath, func(f string, d fs.DirEntry, err error) error {
				if err != nil {
					log.Fatal(err)
				}
				correct, err := filepath.Match("*.fa", f)
				if err != nil {
					return err
				} else if correct {
					qFasta := f

					currLine = lastZ + " " + pairwise + "/" + reference + ".byChrom" + "/" + tFasta + " " + pairwise + "/" + species + ".byChrom" + "/" + qFasta + " --output=" + pairwise + "/" + reference + "/" + species + "." + reference + ".axt --scores=" + matrix + " --format=axt " + par
				}
				return nil
			})

		}
		return nil
	})
	return currLine
}

func usage() {
	fmt.Print(
		"lastZWriter was designed to quickly write out lastZ pairwise inputs by contig where multiple references " +
			"are being used. This function writes a text file where each line is an input for a lastZ pairwise alignment." +
			" It requires that each genome be broken with the 'byname' option of 'faSplit' in kentutils and named in the convention " +
			"'assemblyName.byChrom'. Within the parent directory of the byChrom files, lastZWriter will build a " +
			"directory tree for the outputs of lastZ. At the same level of the byChrom files it will create a set of " +
			"directories with the reference species assembly name. For each species being aligned to that reference " +
			"lastZWriter will specify that lastZ should create an axt file 'aligningSpeciesByChromName.referenceByChromName.axt' " +
			"that refers to the two fasta files used in the alignment. LastZWriter also requires a list of all " +
			"species in the alignment, as well as a separate text file with a list of reference species. " +
			"Matrices are hardcoded absolute paths in this version. Matrices are assigned based on the distance " +
			"between the reference and aligning species from each other as calculated by the PHAST all_dists function, " +
			"which is also a required file. This function can be used directly within the terminal, but would be " +
			"easiest to work with in a shell wrapper where inputs can be referred to in variables. \n" +
			"Usage:\n" +
			"lastZWriter <lastZ install> <path to parent of .byChrom> <speciesList.txt> <referenceList.txt> <allDists.txt> <outFile.txt>" +
			"bedOverlapByWindow takes a sorted bed and counts bp in bed regions within a window size. Default is 5000bp\n" +
			"Usage:\n" +
			"bedOverlapByWindow input.bed chrom.sizes output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

//TODO: write needed species file and reference file
//TODO: write matrix in working dir while this is running and it can be accessed during that run time for the program, leave this as an option where the default is hard-coding a path

func main() {
	var expectedNumArgs int = 7
	var mergeAdjacent *bool = flag.Bool("mergeAdjacent", false, "Merge non-overlapping entries with direct adjacency.")
	var lowMem *bool = flag.Bool("lowMem", false, "Use the low memory algorithm. Requires input file to be pre-sorted.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	bedMerge(infile, outfile, *mergeAdjacent, *lowMem)
}
