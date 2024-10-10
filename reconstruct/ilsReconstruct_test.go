package reconstruct

import (
	"testing"
	"bufio"
	"os"
	"fmt"

	"github.com/vergenlab/gononomics/simulate"
	"github.com/vergenlab/gononomics/exception"
	"github.com/vergenlab/gononomics/fasta/pFasta"
	"github.com/vergenlab/gononomics/fileio"
)

var IlsReconstructTests = []struct {
	PostProbsFileName	string
	ReconsFileName		string
	Precision			float32
	OutFile				string
	Expected			string
}{
	{PostProbsFileName: "testdata/ilsPostProbs.txt",
	ReconsFileName: "testdata/ilsReconsInput.txt",
	Precision: 0.001,
	OutFile: "testdata/ilsRecon.pFa",
	Expected: "testdata/ilsRecon_Expected.pFa"
	TestcasePrecision: 0.001,
	},
}

func TestIlsReconstruct(t *testing.T) {
	var out pFasta.PFasta
	var fileScanner *bufio.Scanner
	var error Error
	var readPostProbs *File
	var readRecons *File
	var postProbFileLines []string
	var reconFileLines []string
	var postProbs [][]float64
	var recons []pFasta.PFasta
	var expected pFasta.PFasta
	

	for _, v := range IlsReconstructTests {

		readPostProbs, err = os.Open(v.PostProbsFileName)
		if err != nil {
			log.Fatalf("%s does not exist.", v.PostProbsFileName)
		}

		fileScanner = bufio.NewScanner(readPostProbs)
		fileScanner.Split(bufio.ScanLines)
		postProbFileLines = []

		for fileScanner.Scan() {
			postProbFileLines = append(postProbFileLines, fileScanner.Text())
		}

		readPostProbs.Close()
		postProbs = []

		for _, filepath := range postProbfileLines {
			postProbs = append(postProbs, wig.Read(filepath))
		}


		readRecons, err = os.Open(v.ReconsFileName)
		if err != nil {
			log.Fatalf("%s does not exist.", v.ReconsFileName)
		}

		fileScanner = bufio.NewScanner(readRecons)
		fileScanner.Split(bufio.ScanLines)
		reconFileLines = []

		for fileScanner.Scan() {
			reconFileLines = append(reconFileLines, fileScanner.Text())
		}

		readRecons.Close()
		recons = []pFasta.PFasta{}

		for _, filepath := range reconFileLines {
			recons = append(recons, pFasta.Read(filepath))
		}


		out = IlsReconstructSeq(PostProbsFileName)
		pFasta.Write(v.OutFile, out)

		expected = pFasta.Read(v.Expected)
		if !pFasta.AllAreEqual(out, expected, v.TestcasePrecision) {
			t.Errorf("Error: ilsReconstruct outfile is not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
	}
}