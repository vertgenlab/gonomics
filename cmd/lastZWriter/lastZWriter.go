package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/lastZWriter"
	"io/fs"
	"log"
	"path/filepath"
)

func MakeArray(lastZ string, pairwise string, speciesListFile string, refListFile string, allDists string, outText string, m bool, mPath string) {
	if !m {
		lastZWriter.buildMatrices(mPath)
	}
	speciesList := fileio.EasyOpen(speciesListFile)
	refList := fileio.EasyOpen(refListFile)
	fileio.EasyCreate(outText)
	var spec, ref string
	var speciesDone, refDone bool
	var parameters []string
	var matrix string
	var dist int
	var allLines []string
	for ref, refDone = fileio.EasyNextRealLine(refList); !refDone; ref, refDone = fileio.EasyNextRealLine(refList) {
		for spec, speciesDone = fileio.EasyNextRealLine(speciesList); !speciesDone; spec, speciesDone = fileio.EasyNextRealLine(speciesList) {
			if spec != ref {
				parameters, matrix, dist = lastZWriter.AlignSetUp(pairwise, spec, ref, allDists, m, mPath)
				allLines = writeFiles(lastZ, pairwise, ref, spec, parameters, matrix, dist, allLines)
			}
		}
	}
	fileio.Write(outText, allLines)
}

func writeFiles(lastZ string, pairwise string, reference string, species string, parameters []string, matrix string, dist int, allLines []string) (lines []string) {
	var currLine string
	par := fmt.Sprintf("%s %s %s %s %s %s %s ", parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], parameters[6])

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
			"lastZWriter will specify that lastZ should create an axt output file 'aligningSpeciesByChromName.referenceByChromName.axt' " +
			"that refers to the two fasta files used in the alignment. LastZWriter also requires a list of all " +
			"species in the alignment, as well as a separate text file with a list of reference species. " +
			"Matrices are hardcoded, absolute paths by default in this version. In the default function matrices are assigned based on the distance " +
			"between the reference and aligning species from each other as calculated by the PHAST all_dists function." +
			"However, the user has the option to specify a bool (option m) as false and provide a path in which they " +
			"would like the needed matrices to be hardcoded. As an alternative to the all_dists function, or if there " +
			"isn't an available tree of the necessary species, this function can also take a file to replace the " +
			"specified allDists file. The first two columns of which would need to be every possible combination of " +
			"your alignment (find an example in testdata directory). This function can be used directly within the " +
			"terminal, but would be easiest to work with in a shell wrapper where inputs can be referred to in variables. \n" +
			"Usage:\n" +
			"lastZWriter [-m=<bool> -mPath=<string>] <lastZ install> <path to parent of .byChrom> <speciesList.txt> <referenceList.txt> <allDists.txt> <outFile.txt>" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 6
	var m *bool = flag.Bool("m", true, "Use existing matrices at hardcoded path.")
	var mPath *string = flag.String("mPath", "", "Path to desired location of created matrices if m = false.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	lastZ := flag.Arg(0)
	pairwiseDir := flag.Arg(1)
	speciesListFile := flag.Arg(2)
	refListFile := flag.Arg(3)
	allDists := flag.Arg(4)
	outText := flag.Arg(5)

	MakeArray(lastZ, pairwiseDir, speciesListFile, refListFile, allDists, outText, *m, *mPath)
}
