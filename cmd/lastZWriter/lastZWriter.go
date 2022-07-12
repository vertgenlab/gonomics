package main

import (
	"encoding/csv"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/lastZWriter"
)

func MakeArray(speciesListFile string, refListFile string, allDists string, outText string, outCsv string) {
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
				parameters, matrix, dist = lastZWriter.AlignSetUp(spec, ref, allDists)
				allLines := writeFiles(ref, spec, parameters, matrix, dist, allLines, csvOut)
			}
		}
	}
	fileio.Write(outText, allLines)

}

func writeFiles(reference string, species string, parameters []string, matrix string, dist int, allLines []string, csvOut *csv.Writer) (lines []string) {
	var currLine string
	lastZ := "/hpc/group/vertgenlab/softwareShared/lastz-master/src/lastz"

	allVars := []string{
		reference, species, string(dist), matrix}

	csvOut.Write(allVars)
	csvOut.Write(parameters)

	currLine = lastZ
	allLines = append(allLines, currLine)

	return allLines

	//TODO:add lines here but write and update files in the wrapper (make array)

}
