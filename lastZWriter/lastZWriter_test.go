package lastZWriter

import (
	"log"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var parClose = []string{"O=600", "E=150", "H=2000", "T=2", "M=254", "K=4500", "L=3000", "Y=15000"}
var parDefault = []string{"O=400", "E=30", "H=2000", "T=1", "M=254", "K=3000", "L=3000", "Y=9400"}
var parFar = []string{"O=400", "E=30", "H=2000", "T=1", "M=50", "K=2200", "L=6000", "Y=3400"}
var mClose = "/hpc/group/vertgenlab/alignmentSupportFiles/human_chimp_v2.mat"
var mDefault = "/hpc/group/vertgenlab/alignmentSupportFiles/default.mat"
var mFar = "/hpc/group/vertgenlab/alignmentSupportFiles/hoxD55.mat"
var mClosePath = "testdata/human_chimp_v2.mat"
var mDefaultPath = "testdata/default.mat"
var mFarPath = "testdata/hoxD55.mat"

var expectedPaths = []string{
	"testdata/refer1.refer2/chr10", "testdata/refer1.refer2/chr11", "testdata/refer1.name1/chr10",
	"testdata/refer1.name1/chr11", "testdata/refer1.name2/chr10", "testdata/refer1.name2/chr11",
	"testdata/refer2.refer1/chr12", "testdata/refer2.refer1/chr13", "testdata/refer2.name1/chr12",
	"testdata/refer2.name1/chr13", "testdata/refer2.name2/chr12", "testdata/refer2.name2/chr13",
}

var parentDirs = []string{
	"testdata/refer1.refer2", "testdata/refer1.name1", "testdata/refer1.name2",
	"testdata/refer2.refer1", "testdata/refer2.name1", "testdata/refer2.name2",
}

var matrices = []string{
	"testdata/default.mat", "testdata/hoxD55.mat", "testdata/human_chimp_v2.mat",
}

func TestAlignSetUp(t *testing.T) {
	pairDir := "testdata"
	m := true
	rList := fileio.Read(pairDir + "/refList.txt")
	specList := fileio.Read(pairDir + "/speciesList.txt")
	for ref := range rList {
		for species := range specList {
			if rList[ref] == specList[species] {
				continue
			}
			par, mat := AlignSetUp(pairDir, specList[species], rList[ref], pairDir+"/allDistsAll.txt", m, "")
			if rList[ref] == "refer1" && specList[species] == "refer2" || rList[ref] == "refer2" && specList[species] == "refer1" {
				for s := range par {
					parMatch := strings.Compare(par[s], parClose[s])
					if parMatch != 0 {
						t.Fatalf("Close parameters mismatch: %v value didn't match", s)
					}
				}
				mMatch := strings.Compare(mat, mClose)
				if mMatch != 0 {
					t.Fatal("Close matrix mismatch")
				}

			} else if rList[ref] == "refer1" && specList[species] == "name1" || rList[ref] == "refer2" && specList[species] == "name1" {
				for s := range par {
					parMatch := strings.Compare(par[s], parDefault[s])
					if parMatch != 0 {
						t.Fatalf("Default parameters mismatch: %v value didn't match", s)
					}
				}
				mMatch := strings.Compare(mat, mDefault)
				if mMatch != 0 {
					t.Fatal("Default matrix mismatch")
				}
			} else if rList[ref] == "refer1" && specList[species] == "name2" || rList[ref] == "refer2" && specList[species] == "name2" {
				for s := range par {
					parMatch := strings.Compare(par[s], parFar[s])
					if parMatch != 0 {
						t.Fatalf("Far parameters mismatch: %v value didn't match", s)
					}
				}
				mMatch := strings.Compare(mat, mFar)
				if mMatch != 0 {
					t.Fatal("Far matrix mismatch")
				}
			} else {
				log.Panicf("Didn't assign any matrix or parameter to this alignment, Ref: %s, Spec: %s", rList[ref], specList[species])
			}
		}
	}
	for p := range expectedPaths {
		_, e := os.Stat(expectedPaths[p])
		if e != nil {
			t.Fatalf("Expected path check returned error %e", e)
		}
	}
	removeThings()
}

func TestAlignSetUp2(t *testing.T) {
	var mPath = "testdata"
	pairDir := "testdata"
	m := false
	rList := fileio.Read(pairDir + "/refList.txt")
	specList := fileio.Read(pairDir + "/speciesList.txt")
	for ref := range rList {
		for species := range specList {
			if rList[ref] == specList[species] {
				continue
			}
			if !m {
				BuildMatrices(mPath)
			}
			par, mat := AlignSetUp(pairDir, specList[species], rList[ref], pairDir+"/allDistsAll.txt", m, mPath)
			if rList[ref] == "refer1" && specList[species] == "refer2" || rList[ref] == "refer2" && specList[species] == "refer1" {
				for s := range par {
					parMatch := strings.Compare(par[s], parClose[s])
					if parMatch != 0 {
						t.Fatalf("Close parameters mismatch: %v value didn't match", s)
					}
				}
				mMatch := strings.Compare(mat, mClosePath)
				if mMatch != 0 {
					t.Fatal("Close matrix mismatch")
				}
			} else if rList[ref] == "refer1" && specList[species] == "name1" || rList[ref] == "refer2" && specList[species] == "name1" {
				for s := range par {
					parMatch := strings.Compare(par[s], parDefault[s])
					if parMatch != 0 {
						t.Fatalf("Default parameters mismatch: %v value didn't match", s)
					}
				}
				mMatch := strings.Compare(mat, mDefaultPath)
				if mMatch != 0 {
					t.Fatal("Default matrix mismatch")
				}
			} else if rList[ref] == "refer1" && specList[species] == "name2" || rList[ref] == "refer2" && specList[species] == "name2" {
				for s := range par {
					parMatch := strings.Compare(par[s], parFar[s])
					if parMatch != 0 {
						t.Fatalf("Far parameters mismatch: %v value didn't match", s)
					}
				}
				mMatch := strings.Compare(mat, mFarPath)
				if mMatch != 0 {
					t.Fatal("Far matrix mismatch")
				}
			} else {
				log.Panicf("Didn't assign any matrix or parameter to this alignment, Ref: %s, Spec: %s", rList[ref], specList[species])
			}
		}
	}
	for p := range expectedPaths {
		_, e := os.Stat(expectedPaths[p])
		if e != nil {
			t.Fatalf("Expected path check returned error %e", e)
		}
	}
	removeThings()
}

func removeThings() {
	for r := range parentDirs {
		os.RemoveAll(parentDirs[r])
	}

	for m := range matrices {
		os.RemoveAll(matrices[m])
	}

}
