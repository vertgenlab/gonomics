package reconstruct

import (
	"bufio"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"os"
	"testing"
	// "github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	// "math/rand"
	// "github.com/vertgenlab/gonomics/dna/pDna"
)

var IlsReconstructTests = []struct {
	PostProbsFileName  string
	ReconsFileName     string
	ChromSizesFileName string
	Precision          float32
	OutFile            string
	Expected           string
	TestcasePrecision  float32
}{
	{PostProbsFileName: "testdata/ilsPostProbs.txt",
		ReconsFileName:     "testdata/ilsReconsInput.txt",
		ChromSizesFileName: "testdata/ilschromSizes.chrom.sizes",
		Precision:          0.001,
		OutFile:            "testdata/ilsRecon.pFa",
		Expected:           "testdata/ilsRecon_Expected.pFa",
		TestcasePrecision:  0.001,
	},
}

func TestIlsReconstruct(t *testing.T) {
	var out []pFasta.PFasta
	var fileScanner *bufio.Scanner
	var err error
	var readPostProbs *os.File
	var readRecons *os.File
	var postProbFileLines []string
	var reconFileLines []string
	var postProbs []map[string]wig.Wig
	var recons []pFasta.PFasta
	var expected []pFasta.PFasta

	for _, v := range IlsReconstructTests {
		readRecons, err = os.Open(v.ReconsFileName)
		if err != nil {
			log.Fatalf("%s does not exist.", v.ReconsFileName)
		}

		fileScanner = bufio.NewScanner(readRecons)
		fileScanner.Split(bufio.ScanLines)
		reconFileLines = make([]string, 0)

		for fileScanner.Scan() {
			reconFileLines = append(reconFileLines, fileScanner.Text())
		}

		readRecons.Close()
		recons = make([]pFasta.PFasta, 0)

		for _, filepath := range reconFileLines {
			recons = append(recons, pFasta.Read(filepath)[0])
		}

		readPostProbs, err = os.Open(v.PostProbsFileName)
		if err != nil {
			log.Fatalf("%s does not exist.", v.PostProbsFileName)
		}

		fileScanner = bufio.NewScanner(readPostProbs)
		fileScanner.Split(bufio.ScanLines)
		postProbFileLines = make([]string, 0)

		for fileScanner.Scan() {
			postProbFileLines = append(postProbFileLines, fileScanner.Text())
		}

		readPostProbs.Close()
		postProbs = make([]map[string]wig.Wig, 0)

		for _, filepath := range postProbFileLines {
			postProbs = append(postProbs, wig.Read(filepath, v.ChromSizesFileName, 0))
		}

		out = []pFasta.PFasta{IlsReconstructSeq(postProbs, recons, v.Precision)}
		pFasta.Write(v.OutFile, out)

		expected = pFasta.Read(v.Expected)
		if !pFasta.AllAreEqual(out, expected, v.TestcasePrecision) {
			t.Errorf("Error: ilsReconstruct outfile is not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
	}
}

// var IdxTests = []struct {
// 	Idx	int64
// } {
// 	{Idx: 4,},
// }

// func TestIdx(t* testing.T) {
// 	for _, v := range IdxTests {
// 		var out pFasta.PFasta
// 		for idx := range v.Idx {
// 			out.Seq = make([]pDna.Float32Base, 0)
// 			out.Name = fmt.Sprintf("chr1")
// 			r := rand.New(rand.NewSource(10+idx))
// 			sum := float32(0)
// 			pos := 0
// 			var aInt float32
// 			var cInt float32
// 			var gInt float32
// 			var tInt float32
// 			var base pDna.Float32Base
// 			for pos < 11 {
// 				aInt = r.Float32()
// 				cInt = r.Float32()
// 				gInt = r.Float32()
// 				tInt = r.Float32()
// 				sum = aInt + cInt + gInt + tInt
// 				fmt.Printf("\nsum %v\n", sum)
// 				base.A = float32(aInt/sum)
// 				base.C = float32(cInt/sum)
// 				base.G = float32(gInt/sum)
// 				base.T = float32(tInt/sum)
// 				pos += 1
// 				out.Seq = append(out.Seq, base)
// 			}
// 			pos = 0
// 			outFileName := fmt.Sprintf("testdata/ilsRecon_%v.pfa", idx)
// 			pFasta.Write(outFileName, []pFasta.PFasta{out})
// 			fmt.Printf("Finish round %v", idx)
// 		}
// 	}

// 	fmt.Print("DONE")
// }
