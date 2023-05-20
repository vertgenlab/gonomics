package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/lastZWriter"
)

func MakeArray(lastZ string, pairwise string, speciesListFile string, refListFile string, allDists string, outText string, m bool, mPath string) {
	if !m {
		lastZWriter.BuildMatrices(mPath)
	}
	speciesList := fileio.Read(speciesListFile)
	refList := fileio.Read(refListFile)
	fileio.EasyCreate(outText)
	var parameters []string
	var matrix string
	var allLines []string
	for ref := range refList {
		for spec := range speciesList {
			match := strings.Compare(speciesList[spec], refList[ref])
			if match != 0 {
				parameters, matrix = lastZWriter.AlignSetUp(pairwise, speciesList[spec], refList[ref], allDists, m, mPath)
				if parameters == nil || matrix == "" {
					log.Fatalf("Reference %s and species %s returned no parameters or matrix.", refList[ref], speciesList[spec])
				}
				allLines = writeFile(lastZ, pairwise, refList[ref], speciesList[spec], parameters, matrix, allLines)
			}
		}
	}
	fileio.Write(outText, allLines)
}

// MakeArraySimple has no allDists, m, or mPath (will not generate matrix), and parameters is directly input as 1 string not generated as []string
func MakeArraySimple(lastZ string, pairwise string, speciesListFile string, refListFile string, parameters string, outText string) {
	speciesList := fileio.Read(speciesListFile)
	refList := fileio.Read(refListFile)
	fileio.EasyCreate(outText)
	var allLines []string
	for ref := range refList {
		for spec := range speciesList {
			match := strings.Compare(speciesList[spec], refList[ref])
			if match != 0 {
				lastZWriter.AlignSetUpSimple(pairwise, speciesList[spec], refList[ref])
				allLines = writeFileSimple(lastZ, pairwise, refList[ref], speciesList[spec], parameters, allLines)
			}
		}
	}
	fileio.Write(outText, allLines)
}

func writeFile(lastZ string, pairwise string, reference string, species string, parameters []string, matrix string, allLines []string) (lines []string) {
	var currLines []string
	par := fmt.Sprintf("%s %s %s %s %s %s %s %s ", parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], parameters[6], parameters[7])

	currLines = fastaFinder(lastZ, pairwise, reference, species, par, matrix)
	allLines = append(allLines, currLines...)
	return allLines
}

// writeFileSimple has no matrix string, and parameters is 1 string not []string
func writeFileSimple(lastZ string, pairwise string, reference string, species string, parameters string, allLines []string) (lines []string) {
	var currLines []string
	currLines = fastaFinderSimple(lastZ, pairwise, reference, species, parameters)
	allLines = append(allLines, currLines...)
	return allLines
}

func fastaFinder(lastZ string, pairwise string, reference string, species string, par string, matrix string) (lines []string) {
	var currLine string
	var theseLines []string
	var tMatches, qMatches, tFiles, qFiles []string
	tPath := filepath.Join(pairwise, reference + ".byChrom")
	qPath := filepath.Join(pairwise, species + ".byChrom")

	if _, e := os.Stat(tPath); os.IsNotExist(e) {
		log.Fatalf("There is no .byChrom directory for the target (reference) species.")
	}
	if _, e := os.Stat(qPath); os.IsNotExist(e) {
		log.Fatalf("There is no .byChrom directory for the query species.")
	}

	tMatches, _ = filepath.Glob(tPath + "/*.fa")
	qMatches, _ = filepath.Glob(qPath + "/*.fa")

	for tF := range tMatches {
		_, tName := filepath.Split(tMatches[tF])
		tFiles = append(tFiles, tName)
	}
	for qF := range qMatches {
		_, qName := filepath.Split(qMatches[qF])
		qFiles = append(qFiles, qName)
	}
	for t := range tFiles {
		tName := strings.TrimSuffix(tFiles[t], ".fa")
		for q := range qFiles {
			qName := strings.TrimSuffix(qFiles[q], ".fa")
			currLine := lastZ + " " + filepath.Join(
				pairwise,
				reference + ".byChrom",
				tFiles[t]) + " " + filepath.Join(
				pairwise,
				species + ".byChrom",
				qFiles[q]) + " --output=" + filepath.Join(
				pairwise,
				reference + "." + species,
				tName,
				qName + "." + tName + ".axt") + " --scores=" + matrix + " --action:target=multiple" + " --format=axt " + par
			theseLines = append(theseLines, currLine)
		}
	}

	if theseLines == nil {
		log.Fatal("No lines to write to file")
	}

	return theseLines
}

// fastaFinderSimple has no matrix string
// fastaFinderSimple has a different output file name structure: ref.species/qName/tName.qName.axt
func fastaFinderSimple(lastZ string, pairwise string, reference string, species string, par string) (lines []string) {
	var currLine string
	var theseLines []string
	var tMatches, qMatches, tFiles, qFiles []string
	tPath := filepath.Join(pairwise, reference + ".byChrom")
	qPath := filepath.Join(pairwise, species + ".byChrom")

	if _, e := os.Stat(tPath); os.IsNotExist(e) {
		log.Fatalf("There is no .byChrom directory for the target (reference) species.")
	}
	if _, e := os.Stat(qPath); os.IsNotExist(e) {
		log.Fatalf("Error: There is no .byChrom directory for the query species.")
	}

	tMatches, _ = filepath.Glob(tPath + "/*.fa")
	qMatches, _ = filepath.Glob(qPath + "/*.fa")

	for tF := range tMatches {
		_, tName := filepath.Split(tMatches[tF])
		tFiles = append(tFiles, tName)
	}
	for qF := range qMatches {
		_, qName := filepath.Split(qMatches[qF])
		qFiles = append(qFiles, qName)
	}
	for t := range tFiles {
		if strings.HasSuffix(tFiles[t], ".fa") {
			tName := strings.TrimSuffix(tFiles[t], ".fa")
		} else if strings.HasSuffix(tFiles[t], ".fasta" {
			tName := strings.TrimSuffix(tFiles[t], ".fasta")
		} else {
			log.Errorf("Error: there is a bug in the glob command, we found files of other suffixes!!")
		}
		for q := range qFiles {
			qName := strings.TrimSuffix(qFiles[q], ".fa")
			// currLine reflects that fastaFinderSimple has no matrix string
			// currLine also reflects that fastaFinderSimple has a different output file name structure: ref.species/qName/tName.qName.axt
			currLine := lastZ + " " + filepath.Join(
				pairwise,
				reference + ".byChrom",
				tFiles[t]) + " " + filepath.Join(
				pairwise,
				species + ".byChrom",
				qFiles[q]) + " --output=" + filepath.Join(
				pairwise,
				reference + "." + species,
				qName,
				tName + "." + qName + ".axt") + " --action:target=multiple" + " --format=axt " + par
			theseLines = append(theseLines, currLine)
		}
	}

	if theseLines == nil {
		log.Fatal("Error: No lines to write to file")
	}

	return theseLines
}

func usage() {
	fmt.Print(
		"lastZWriter was designed to write out lastZ pairwise inputs by contig where multiple references " +
			"are being used. This function writes a text file where each line is an input for a lastZ pairwise alignment." +
			" It requires that each genome be broken with the 'byname' option of 'faSplit' in kentutils and named in the convention " +
			"'assemblyName.byChrom'. Within the parent directory of the byChrom directories, lastZWriter will build a " +
			"directory tree for the outputs of lastZ. At the same level of the byChrom files it will create a set of " +
			"directories with the naming convention reference.aligned within which will be directories for each contig " +
			"of the reference genome. For each species being aligned to that reference lastZWriter will specify that " +
			"lastZ should create an axt output file 'pairwiseDir/ref.species/referenceContig/aligningSpeciesByChromName.referenceByChromName.axt' " +
			"that refers to the two fasta files used in the alignment. LastZWriter also requires a list of all " +
			"species in the alignment, as well as a separate text file with a list of reference species. " +
			"Matrices are hardcoded absolute paths by default in this version. In the default function matrices are assigned based on the distance " +
			"between the reference and aligning species from each other as calculated by the PHAST all_dists function." +
			"However, the user has the option to specify a bool (option m) as false and provide a path in which they " +
			"would like the needed matrices to be hardcoded. As an alternative to the all_dists function, or if there " +
			"isn't an available tree of the necessary species, this function can also take a file to replace the " +
			"specified allDists file. The first two columns of which would need to be every possible combination of " +
			"their alignment (find an example in gonomics/lastZWriter/testdata directory). This function can be used directly within the " +
			"terminal, but would be easiest to work with in a shell wrapper where inputs can be referred to in variables. \n" +
			"\n" +
			"Instead of regular lastZWriter, the user has the option to run simple lastZWriter, " +
			"which accepts a parameter string directly and does not generate matrices. " +
			"One usage simple lastZWriter is appropriate for is aligning relatively short sequences with whole chromosomes from the same species. " +
			"To run simple lastZWriter, set option simple to true, and provide the parameter string in the option parameters if needed. " +
			"For simple lastZWriter, allDists is still a required input but will not be used, so it can be a dummy input, e.g. 'allDists.txt'. \n" +
			"\n" +
			"Usage:\n" +
			"lastZWriter [-m=<bool> -mPath=<string> -simple=<bool> -parameters=<string>] <lastZ install> <path to parent of .byChrom> <speciesList.txt> <referenceList.txt> <allDists.txt> <outFile.txt> \n" +
			"options:\n")
	flag.PrintDefaults()
}

//TODO: edit usage after editing out make query subdir function

func main() {
	var expectedNumArgs int = 6
	var m *bool = flag.Bool("m", true, "Use existing matrices at hardcoded path.")
	var mPath *string = flag.String("mPath", "", "Path to desired location of created matrices if m = false.")
	// to run simple lastZWriter, need 2 options (simple *bool, parameters *string). Do not combine the 2 options because want empty ("") parameters to work too
	var simple *bool = flag.Bool("simple", false, "Instead of regular lastZWriter, run simple lastZWriter.")
	var parameters *string = flag.String("parameters", "", "Instead of regular lastZWriter, run simple lastZWriter with the specified parameter string.")

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

	if *simple {
		MakeArraySimple(lastZ, pairwiseDir, speciesListFile, refListFile, *parameters, outText)
	} else {
		MakeArray(lastZ, pairwiseDir, speciesListFile, refListFile, allDists, outText, *m, *mPath)
	}
}
