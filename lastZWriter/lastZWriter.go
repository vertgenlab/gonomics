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
//func makeArray(speciesListFile string, refListFile string, allDists string) (par []string, mat string, dis int) { //returns for each ref and spec pairing that it is given, the targetDir, the queryDir, t, q, outDir, tName, qName, matrix, paramters
//	speciesList := fileio.EasyOpen(speciesListFile)
//	refList := fileio.EasyOpen(refListFile)
//	var spec, ref string
//	var speciesDone, refDone bool
//	var parameters []string
//	var matrix string
//	var dist int
//	for ref, refDone = fileio.EasyNextRealLine(refList); !refDone; ref, refDone = fileio.EasyNextRealLine(refList) {
//		for spec, speciesDone = fileio.EasyNextRealLine(speciesList); !speciesDone; spec, speciesDone = fileio.EasyNextRealLine(speciesList) {
//			if spec != ref {
//				parameters, matrix, dist = AlignSetUp(spec, ref, allDists)
//				return parameters, matrix, dist
//			}
//		}
//	}
//	return
//}

//AlignSetUp takes in a the path to the parent directory where all the individual fastas for alignment are being held.
//In this context the directory below "pairwise" is labelled by a species name, which then contains all the fastas for
//alignment. It also takes a single aligning species and reference species, as well as a text file that describes the
//distance between all species in the alignment from all other species in the alignment. This file is make with
//Phylogenetic Analysis with Space/Time Models or PHAST all_dists function. AlignSetUp then calls its helper functions
//and returns the results of findParameters.
func AlignSetUp(pairwise string, species string, reference string, allDists string, m bool, mPath string) (par []string, mat string, dis float64) {
	outDir := pairwise + "/" + reference + "." + species
	makeOutDir(pairwise, outDir, reference, species)
	parameters, matrix, dist := findParameters(reference, species, allDists, m, mPath)
	return parameters, matrix, dist
}

//makeOutDir creates the file directory tree where the output of all of the alignments will go by first creating
//the directory labelled with the name of the reference and the species being aligned. It then passes off to a
//helper function makeTargetSubDir
func makeOutDir(pairwise string, outDir string, r string, s string) {
	tDir := pairwise + "/" + r + ".byChrom"
	if _, e := os.Stat(outDir); os.IsNotExist(e) {
		err := os.Mkdir(outDir, 666)
		if err != nil {
			log.Fatal(err)
		}
	}
	makeTargetSubDir(tDir, outDir, pairwise, s)
}

//makeTargetSubDir creates the next directory layer below makeOutDir which contains all alingments to a
//single reference against any other species.
func makeTargetSubDir(path string, outDir string, pairwise string, s string) {
	var tr, trName string
	//qDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + s + ".byChrom"
	qDir := pairwise + "/" + s + ".byChrom"
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
func findParameters(reference string, species string, distsFile string, m bool, mPath string) (par []string, matrix string, dis float64) {
	var words []string
	var answer []string
	var dist float64
	var mat string
	dists := fileio.EasyOpen(distsFile)
	for line, done := fileio.EasyNextRealLine(dists); !done; line, done = fileio.EasyNextRealLine(dists) {
		words = strings.Split(line, " ")
		if words[0] == reference && words[1] == species {
			if words[2] == "close" {
				answer = append(answer, "O=600", "E=150", "T=2", "M=254", "K=4500", "L=3000", "Y=15000")
				if m {
					matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/human_chimp_v2.mat"
				} else {
					matrix = mPath + "/human_chimp_v2.mat"
				}
				dist = 1
			} else if words[2] == "far" {
				answer = append(answer, "O=400", "E=30", "T=1", "M=50", "K=2200", "L=6000", "Y=3400")
				matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/hoxD55.mat"
				dist = 3
			} else if words[2] == "default" {
				answer = append(answer, "O=400", "E=30", "T=1", "M=254", "K=3000", "L=3000", "Y=9400")
				matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/default.mat"
				dist = 2
			} else {
				dist = common.StringToFloat64(words[2])
				//dist = d
				switch {
				case dist <= 0.2: //closest
					answer = append(answer, "O=600", "E=150", "T=2", "M=254", "K=4500", "L=3000", "Y=15000")
					matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/human_chimp_v2.mat"
				case dist >= 0.7: //farthest
					answer = append(answer, "O=400", "E=30", "T=1", "M=50", "K=2200", "L=6000", "Y=3400")
					matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/hoxD55.mat"
				default: //executive decision to set M to 254
					answer = append(answer, "O=400", "E=30", "T=1", "M=254", "K=3000", "L=3000", "Y=9400")
					matrix = "/hpc/group/vertgenlab/alignmentSupportFiles/default.mat"
				}
			}
		}
	}

	return answer, mat, dist
}

func buildMatrices(mPath string) {
	err := os.Mkdir(mPath, 666)
	if err != nil {
		log.Panic(err)
	}
	closeRec := "A\tC\tG\tT\n" +
		"A\t90\t-330\t-236\t-356\n" +
		"C\t-330\t100\t-318\t-236\n" +
		"G\t-236\t-318\t100\t-330\n" +
		"T\t-356\t-236\t-330\t90"
	fileio.Write(mPath+"/human_chimp_v2.mat", closeRec)

	defaultRec := "A\tC\tG\tT\n" +
		"A\t91\t-114\t-31\t-123\n" +
		"C\t-114\t100\t-125\t-31\n" +
		"G\t-31\t-125\t100\t-114\n" +
		"T\t-123\t-31\t-114\t91"
	fileio.Write(mPath+"/default.mat", defaultRec)

	farRec := "A\tC\tG\tT\n" +
		"A\t91\t-90\t-25\t-100\n" +
		"C\t-90\t100\t-100\t-25\n" +
		"G\t-25\t-100\t100\t-90\n" +
		"T\t-100\t-25\t-90\t91"
	fileio.Write(mPath+"/hoxD55.mat", farRec)
}
