package reconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	//DEBUG: "fmt"
)

//PrimateRecon returns a human-biased conservative human-chimp ancestor estimate from an input multiFa alignment in the order Human-Chimp-Bonoboo-Orangutan-Gorilla. MessyToN, when true, turns ambiguous bases to Ns in the output.
func PrimateRecon(records []fasta.Fasta, messyToN bool) fasta.Fasta {
	answer := fasta.Fasta{Name: "Human_Chimp_Ancestor"}

	var outputBase dna.Base

	if len(records) != 5 {
		log.Fatalf("Wrong number of sequences, expecting five, found %d.\n", len(records))
	}

	//confirm alignment lengths are all the same
	firstLength := len(records[0].Seq)
	for i := 1; i < len(records); i++ {
		if len(records[i].Seq) != firstLength {
			log.Fatalf("Sequence %d is the wrong length.\n", i+1)
		}
	}

	var Ncount int = 0
	var allMatch int = 0
	var humanMatchChimporBonobo int = 0
	var humanChange int = 0
	var gorillaVote int = 0
	var orangutanVote int = 0
	var messyBase int = 0

	human := records[0]
	bonobo := records[1]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	for j := 0; j < len(human.Seq); j++ {
		outputBase = human.Seq[j]
		if human.Seq[j] == dna.N { //output seq is N if human is N
			Ncount++
			outputBase = human.Seq[j]
		} else if humanInsertion(records, j) { //output seq is gap if human is insertion(see helper function)
			outputBase = dna.Gap
		} else if human.Seq[j] != dna.Gap && (chimp.Seq[j] == dna.Gap && bonobo.Seq[j] == dna.Gap) { //if there is sequence in humans and either gorilla or orangutan, but no sequence in chimp and bonobo, we output N, because we don't want to estimate the longer human to gorilla branch.
			if messyToN {
				outputBase = dna.N
			} else {
				outputBase = human.Seq[j]
			}
		} else if gorilla.Seq[j] == dna.Gap && orangutan.Seq[j] == dna.Gap { //no sequence outside the HCA (no gorilla or orangutan, ancestral state can't be determined)
			if messyToN {
				outputBase = dna.N
			} else {
				outputBase = human.Seq[j]
			}
		} else if human.Seq[j] == chimp.Seq[j] && human.Seq[j] == bonobo.Seq[j] { //if human matches chimp and bonobo, output is human
			outputBase = human.Seq[j]
			allMatch++
		} else if (human.Seq[j] == chimp.Seq[j] || human.Seq[j] == bonobo.Seq[j]) && human.Seq[j] != dna.Gap { //if human is a base and matches chimp or bonobo, output is human
			outputBase = human.Seq[j]
			humanMatchChimporBonobo++
		} else if (chimp.Seq[j] == bonobo.Seq[j] && (chimp.Seq[j] == gorilla.Seq[j] || chimp.Seq[j] == orangutan.Seq[j])) && (chimp.Seq[j] != dna.N && chimp.Seq[j] != dna.Gap) { //human is different from ancestor, chimp and bonobo agree with ancestor.
			outputBase = chimp.Seq[j]
			humanChange++
		} else if (human.Seq[j] == gorilla.Seq[j] || chimp.Seq[j] == gorilla.Seq[j] || bonobo.Seq[j] == gorilla.Seq[j]) && (gorilla.Seq[j] != dna.N && gorilla.Seq[j] != dna.Gap) { //more disagreement, but gorilla defines the ancestor in this case
			outputBase = gorilla.Seq[j]
			gorillaVote++
		} else if (human.Seq[j] == orangutan.Seq[j] || chimp.Seq[j] == orangutan.Seq[j] || bonobo.Seq[j] == orangutan.Seq[j] || gorilla.Seq[j] == orangutan.Seq[j]) && (orangutan.Seq[j] != dna.N && orangutan.Seq[j] != dna.Gap) {
			outputBase = orangutan.Seq[j]
			orangutanVote++
		} else {
			if human.Seq[j] != dna.Gap && !messyToN {
				outputBase = human.Seq[j]
			} else {
				outputBase = dna.N
			}
			messyBase++
		}
		answer.Seq = append(answer.Seq, outputBase)
	}
	return answer
}

func humanInsertion(records []fasta.Fasta, j int) bool {
	//true if the human sequence is the only one present at a position
	human := records[0]
	bonobo := records[1]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	if (human.Seq[j] != dna.Gap) && (chimp.Seq[j] == dna.Gap && bonobo.Seq[j] == dna.Gap && gorilla.Seq[j] == dna.Gap && orangutan.Seq[j] == dna.Gap) {
		return true
	}
	return false
}
