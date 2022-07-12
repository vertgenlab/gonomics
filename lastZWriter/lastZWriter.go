package lastZWriter

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io/fs"
	"log"
	"os"
	"path/filepath"
	"strings"
)

//TODO: makeArray should be the cmd, needs to go back and find the target and query dirs, recover tName and qName
//TODO: makeArray needs to use the findParameters return (AlignSetUp return) to make the tsv with a helper in the cmd file, rather than just returning them. This function doesnt need a return.

func makeArray(speciesListFile string, refListFile string, allDists string) (par []string, mat string, dis int) { //returns for each ref and spec pairing that it is given, the targetDir, the queryDir, t, q, outDir, tName, qName, matrix, paramters
	speciesList := fileio.EasyOpen(speciesListFile)
	refList := fileio.EasyOpen(refListFile)
	var spec, ref string
	var speciesDone, refDone bool
	var parameters []string
	var matrix string
	var dist int
	for ref, refDone = fileio.EasyNextRealLine(refList); !refDone; ref, refDone = fileio.EasyNextRealLine(refList) {
		for spec, speciesDone = fileio.EasyNextRealLine(speciesList); !speciesDone; spec, speciesDone = fileio.EasyNextRealLine(speciesList) {
			if spec != ref {
				parameters, matrix, dist = AlignSetUp(spec, ref, allDists)
				return parameters, matrix, dist
			}
		}
	}
	return
}

//AlignSetUp takes in a single aligning species and reference species, as well as a text file that describes the
//distance between all species in the alignment from all other species in the alignment. This file is make with
//Phylogenetic Analysis with Space/Time Models or PHAST all_dists function. AlignSetUp then calls its helper functions
//and returns the results of findParameters.
func AlignSetUp(species string, reference string, allDists string) (par []string, mat string, dis int) {
	outDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + reference + "." + species
	makeOutDir(outDir, reference, species)
	parameters, matrix, dist := findParameters(reference, species, allDists)
	return parameters, matrix, dist
}

//makeOutDir creates the file directory tree where the output of all of the alignments will go by first creating
//the directory labelled with the name of the reference and the species being aligned. It then passes off to a
//helper function makeTargetSubDir
func makeOutDir(outDir string, r string, s string) {
	tDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + r + ".byChrom"
	if _, e := os.Stat(outDir); os.IsNotExist(e) {
		err := os.Mkdir(outDir, 666)
		if err != nil {
			log.Fatal(err)
		}
	}
	makeTargetSubDir(tDir, outDir, s)
}

//makeTargetSubDir creates the next directory layer below makeOutDir which contains all alingments to a
//single reference against any other species.
func makeTargetSubDir(path string, outDir string, s string) {
	var tr, trName string
	qDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + s + ".byChrom"
	filepath.WalkDir(path, func(f string, d fs.DirEntry, err error) error {
		if err != nil {
			log.Fatal(err)
		}
		matched, err := filepath.Match("*.fa", f)
		if err != nil {
			return err
		} else if matched {
			tr = f
			trName = strings.TrimSuffix(tr, ".fa")
			err1 := os.Mkdir(outDir+"/"+trName, 666)
			if err1 != nil {
				log.Fatal(err1)
			}
			parentDir := outDir + "/" + trName
			makeQuerySubDir(qDir, parentDir)
		}
		return nil
	})
}

//makeQuerySubDir makes the second layer of subdirectories within both outDir and targetDir which holds the
//alignment output of a reference (it's parent directory name) against a species (the directory name).
//These directories will remain empty until the array created by these functions is run.
func makeQuerySubDir(path string, pDir string) {
	var qu, quName string
	filepath.WalkDir(path, func(f string, d fs.DirEntry, err error) error {
		if err != nil {
			log.Fatal(err)
		}
		matched, err := filepath.Match("*.fa", f)
		if err != nil {
			return err
		} else if matched {
			qu = f
			quName = strings.TrimSuffix(qu, ".fa")
			err := os.Mkdir(pDir+"/"+quName, 666)
			if err != nil {
				log.Fatal(err)
			}
		}
		return nil
	})
}

//findParameters takes in one pair of species to be alignment (one of them the ref and the other not) and determines
//parameter values for lastZ alignment based on the distance between the species as determined by phast allDists.
//a species can be classified into three categories, closest, where the matrix and parameters will match those used
//for a chimp to human alignment with lastZ, middle, where the parameters will be default for lastZ and the matrix and
//far, where the parameters will be set to the most distant alignment parameters and the matrix will be set to HoxD55.
func findParameters(reference string, species string, distsFile string) (par []string, matrix string, dis int) {
	var words []string
	var answer []string
	var dist int
	var mat string
	dists := fileio.EasyOpen(distsFile)
	for line, done := fileio.EasyNextRealLine(dists); !done; line, done = fileio.EasyNextRealLine(dists) {
		words = strings.Split(line, " ")
		//TODO: is it more efficient to do if words[0] == reference and then nest another if of words[1] == species?
		if words[0] == reference && words[1] == species {
			d := common.StringToFloat64(words[2])
			dist = d
			switch {
			case d <= 0.2: //closest
				answer = append(answer, "600", "150", "2", "254", "4500", "3000", "15000")
				matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/human_chimp_v2.mat"
			case d >= 0.7: //farthest
				answer = append(answer, "400", "30", "1", "50", "2200", "6000", "3400")
				matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/hoxD55.mat"
			default: //executive decision to set M to 254
				answer = append(answer, "400", "30", "1", "254", "3000", "3000", "9400")
				matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/default.mat"
			} //TODO: hard code matrices and give an option of reading in a matrix
		}
	}
	return answer, mat, dist
}
