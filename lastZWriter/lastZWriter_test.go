package lastZWriter

import (
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"testing"
)

var parClose = []string{"O=600", "E=150", "T=2", "M=254", "K=4500", "L=3000", "Y=15000"}
var parDefault = []string{"O=400", "E=30", "T=1", "M=254", "K=3000", "L=3000", "Y=9400"}
var parFar = []string{"O=400", "E=30", "T=1", "M=50", "K=2200", "L=6000", "Y=3400"}
var mClose = "/hpc/group/vertgenlab/alignmentSupportFiles/human_chimp_v2.mat"
var mDefault = "/hpc/group/vertgenlab/alignmentSupportFiles/default.mat"
var mFar = "/hpc/group/vertgenlab/alignmentSupportFiles/hoxD55.mat"

func TestAlignSetUp(t *testing.T) {
	pairDir := "testdata"
	m := true
	rList := fileio.EasyOpen(pairDir + "/refList.txt")
	specList := fileio.EasyOpen(pairDir + "/speciesList.txt")
	for ref, rDone := fileio.EasyNextRealLine(rList); !rDone; ref, rDone = fileio.EasyNextRealLine(specList) {
		for species, sDone := fileio.EasyNextRealLine(specList); !sDone; species, sDone = fileio.EasyNextRealLine(specList) {
			if ref == species {
				continue
			}
			par, mat := AlignSetUp(pairDir, species, ref, pairDir+"/allDistsAll.txt", m, "")

			//log.Print(par)
			//log.Print(mat)
			//log.Print(ref)
			//log.Print(species)

			if ref == "refer1" && species == "refer2" || ref == "refer2" && species == "refer1" {
				for s := range par {
					parMatch := strings.Compare(par[s], parClose[s])

					log.Print(s)
					log.Print(par[s])
					log.Print(parClose[s])

					if parMatch != 0 {
						t.Errorf("Close parameters mismatch: %v value didn't match", s)
					}
					mMatch := strings.Compare(mat, mClose)
					if mMatch != 0 {
						t.Error("Close matrix mismatch")
					}
				}
			} else if ref == "refer1" && species == "name1" || ref == "refer2" && species == "name1" {
				for s := range par {
					parMatch := strings.Compare(par[s], parDefault[s])

					log.Print(s)
					log.Print(par[s])
					log.Print(parDefault[s])

					if parMatch != 0 {
						t.Errorf("Defualt parameters mismatch: %v value didn't match", s)
					}
					mMatch := strings.Compare(mat, mDefault)
					if mMatch != 0 {
						t.Error("Defualt matrix mismatch")
					}
				}
			} else if ref == "refer1" && species == "name2" || ref == "refer2" && species == "name2" {
				for s := range par {
					parMatch := strings.Compare(par[s], parFar[s])

					log.Print(s)
					log.Print(par[s])
					log.Print(parFar[s])

					if parMatch != 0 {
						t.Errorf("Far parameters mismatch: %v value didn't match", s)
					}
					mMatch := strings.Compare(mat, mFar)
					if mMatch != 0 {
						t.Error("Far matrix mismatch")
					}
				}
			} else {
				log.Print("Didn't enter any if/elses")
			}
		}
	}
}

//TODO: do an os.IsExist for the paths that the directory tree should have
