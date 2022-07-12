package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/lastZWriter"
)

func main(speciesListFile string, refListFile string, allDists string, outText string, outTsv string) {
	speciesList := fileio.EasyOpen(speciesListFile)
	refList := fileio.EasyOpen(refListFile)
	var spec, ref string
	var speciesDone, refDone bool
	var parameters []string
	var matrix string
	for ref, refDone = fileio.EasyNextRealLine(refList); !refDone; ref, refDone = fileio.EasyNextRealLine(refList) {
		for spec, speciesDone = fileio.EasyNextRealLine(speciesList); !speciesDone; spec, speciesDone = fileio.EasyNextRealLine(speciesList) {
			if spec != ref {
				parameters, matrix = lastZWriter.AlignSetUp(spec, ref, allDists)
				writeFiles(ref, spec, parameters, matrix)
			}
		}
	}
}

func writeFiles(referece string, species string, parameters []string, matrix string, outText string, outTsv string) {
	//TODO: write TSV with directories involved?
	//TODO: write into text file

}
