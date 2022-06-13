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

func makeArray(speciesListFile string, refListFile string, allDists string) (par []string, mat string) { //returns for each ref and spec pairing that it is given, the targetDir, the queryDir, t, q, outDir, tName, qName, matrix, paramters
	speciesList := fileio.EasyOpen(speciesListFile)
	refList := fileio.EasyOpen(refListFile)
	var spec, ref string
	var speciesDone, refDone bool
	var parameters []string
	var matrix string
	for ref, refDone = fileio.EasyNextRealLine(refList); !refDone; ref, refDone = fileio.EasyNextRealLine(refList) {
		for spec, speciesDone = fileio.EasyNextRealLine(speciesList); !speciesDone; spec, speciesDone = fileio.EasyNextRealLine(speciesList) {
			if spec != ref {
				parameters, matrix = AlignSetUp(spec, ref, allDists)
				return parameters, matrix
			}
		}
	}
	return
}

func AlignSetUp(species string, reference string, allDists string) (par []string, mat string) {
	outDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + reference + "." + species
	output, t, q, tName, qName := makeOutDir(outDir, reference, species)
	parameters, matrix := findParameters(reference, species, allDists)
	return parameters, matrix
}

//findParameters takes in one pair of species to be alignment (one of them the ref and the other not) and determines
//parameter values for lastZ alignment based on the distance between the species as determined by phast allDists.
//a species can be classified into three categories, closest, where the matrix and parameters will match those used
//for a chimp to human alignment with lastZ, middle, where the parameters will be default for lastZ and the matrix and
//far, where the parameters will be set to the most distant alignment parameters and the matrix will be set to HoxD55.
func findParameters(reference string, species string, distsFile string) (par []string, matrix string) {
	var words []string
	var answer []string
	var mat string
	dists := fileio.EasyOpen(distsFile)
	for line, done := fileio.EasyNextRealLine(dists); !done; line, done = fileio.EasyNextRealLine(dists) {
		words = strings.Split(line, " ")
		//TODO: is it more efficient to do if words[0] == reference and then nest another if of words[1] == species?
		if words[0] == reference && words[1] == species {
			d := common.StringToFloat64(words[2])
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
	return answer, mat
}

func makeOutDir(outDir string, r string, s string) (out string, target string, query string, tName string, qName string) {
	tDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + r + ".byChrom"
	//TODO: make a variable that stores the majority of the path
	err := os.Chdir("/hpc/group/vertgenlab/vertebrateConservation/pairwise/")
	if err != nil {
		log.Fatal(err)
	}
	if _, e := os.Stat(outDir); os.IsNotExist(e) {
		err2 := os.Mkdir(outDir, 666)
		if err2 != nil {
			log.Fatal(err2)
		}
	}
	tr, qu, trName, quName := makeTargetSubDir(tDir, outDir, s)

	return out, tr, qu, trName, quName
}
func makeTargetSubDir(path string, outDir string, s string) (target string, query string, tName string, qName string) {
	var tr, qu, trName, quName string
	qDir := "/hpc/group/vertgenlab/vertebrateConservation/pairwise/" + s + ".byChrom"
	err := os.Chdir(outDir)
	if err != nil {
		log.Fatal(err)
	}
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
			err := os.Mkdir(tName, 666)
			if err != nil {
				log.Fatal(err)
			}
			qu, quName = makeQuerySubDir(qDir, tName, outDir)
			//TODO: return all returns for the larger function, may have to make a return section in filepath.WalkDir func
		}
		return nil
	})
	return tr, qu, trName, quName
}

func makeQuerySubDir(path string, tName string, outDir string) (q string, qName string) {
	var qu, quName string
	err := os.Chdir(outDir + "/" + tName)
	if err != nil {
		log.Fatal(err)
	}
	filepath.WalkDir(path, func(f string, d fs.DirEntry, err error) error {
		if err != nil {
			log.Fatal(err)
		}
		matched, err := filepath.Match("*.fa", f)
		if err != nil {
			return err
		} else if matched {
			qu = f
			quName = strings.TrimSuffix(q, ".fa")
			err := os.Mkdir(qName, 666)
			if err != nil {
				log.Fatal(err)
			}
		}
		return nil
	})

	return qu, quName
}

//TODO:I could return nothing and just create the directory tree and in the command I could access them and build each line that way
