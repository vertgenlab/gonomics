package motif

/*
import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)


func Match(pm PositionMatrix, record fasta.Fasta, chromName string, cutoff float64, outfile string) {
	var refPos int
	var err error
	out := fileio.EasyCreate(outfile)


	for alnPos := 0; alnPos < len(record.Seq); alnPos++ {


		if record.Seq[alnPos] != dna.Gap {
			refPos++
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// reuse var output for memory
func scoreWindowPosStrand(pm PositionMatrix, seq []dna.Base, alnStart int, chromName string, output bed.Bed) {
	var currAlnPos = alnStart
	var answer float64
	var motifPos int = 0
	for motifPos < len(pm.Mat[0]) {
		switch seq[currAlnPos] {
		case dna.Gap:
			//we advance in alnPos at the end of the loop to check the next position but do not advance the motifPosition.
		case dna.A:
			answer += pm.Mat[0][motifPos]
			motifPos++
		case dna.C:
			answer += pm.Mat[1][motifPos]
			motifPos++
		case dna.G:
			answer += pm.Mat[2][motifPos]
			motifPos++
		case dna.T:
			answer += pm.Mat[3][motifPos]
			motifPos++
		case dna.N:
			return //this will just ignore this position of the sequence and move on.
		default:
			log.Fatalf("Unrecognized base. Cannot score window.")
			return
		}
		currAlnPos++
	}

}
*/
