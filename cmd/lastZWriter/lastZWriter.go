package main

import (
	"encoding/csv"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/lastZWriter"
	"io/fs"
	"log"
	"path/filepath"
)

//TODO: hardcode matrices in the lastZWriter package
func MakeArray(pairwise string, speciesListFile string, refListFile string, allDists string, outText string, outCsv string) {
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
				allLines = writeFiles(pairwise, ref, spec, parameters, matrix, dist, allLines, csvOut)
			}
		}
	}
	fileio.Write(outText, allLines)

}

func writeFiles(pairwise string, reference string, species string, parameters []string, matrix string, dist int, allLines []string, csvOut *csv.Writer) (lines []string) {
	var currLine string
	par := fmt.Sprintf("%s %s %s %s %s %s %s ", parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], parameters[6])
	allVars := []string{
		reference, species, string(dist), matrix, par}
	csvOut.Write(allVars)

	currLine = fastaFinder(pairwise, reference, species, par, matrix)
	allLines = append(allLines, currLine)

	return allLines
}

func fastaFinder(pairwise string, reference string, species string, par string, matrix string) (line string) {
	var currLine string
	tPath := pairwise + "/" + reference
	qPath := pairwise + "/" + species
	lastZ := "/hpc/group/vertgenlab/softwareShared/lastz-master/src/lastz"

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

func main() {

}
