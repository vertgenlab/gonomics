package reconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

//ChimpRecon is like QuickPrimateRecon but returns a chimp-biased ancestor estimate from an input multiFa in the order Human-Bonobo-Chimp-Orangutan-Gorilla.
//This is the same alignment used for the human-biased HCA estimate. However, bonobo is not used in this estimation.
func ChimpAncestorRecon(records []fasta.Fasta, messyToN bool) fasta.Fasta {
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

	//confirm record names
	if records[0].Name != "hg38" {
		log.Fatalf("First record in multiFa must be named hg38. Found %s.", records[0].Name)
	}
	if records[1].Name != "panPan2" {
		log.Fatalf("Second record in multiFa must be named panPan2. Found %s.", records[1].Name)
	}
	if records[2].Name != "panTro6" {
		log.Fatalf("Third record in multiFa must be named panTro6. Found %s.", records[2].Name)
	}
	if records[3].Name != "gorGor5" {
		log.Fatalf("Fourth record in multiFa must be named gorGor5. Found %s.", records[3].Name)
	}
	if records[4].Name != "ponAbe3" {
		log.Fatalf("Fifth record in multiFa must be named ponAbe3. Found %s.", records[4].Name)
	}

	var Ncount int = 0
	var allMatch int = 0
	var chimpChange int = 0
	var gorillaVote int = 0
	var orangutanVote int = 0
	var messyBase int = 0

	human := records[0]
	//bonobo := records[1]//This is unused for the chimp estimation.
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	for j := 0; j < len(human.Seq); j++ {
		if chimp.Seq[j] == dna.N {
			Ncount++
			outputBase = chimp.Seq[j]
		} else if chimpInsertion(records, j) { //output base is gap if the chimp sequence is an insertion.
			outputBase = dna.Gap
		} else if chimp.Seq[j] != dna.Gap && human.Seq[j] == dna.Gap { //if there is sequence in chimp and either gorilla or orangutan, but no sequence in human, we output N, because we don't want to estimate the longer human to gorilla branch.
			if messyToN {
				outputBase = dna.N
			} else {
				outputBase = chimp.Seq[j]
			}
		} else if gorilla.Seq[j] == dna.Gap && orangutan.Seq[j] == dna.Gap { //no sequence outside the HCA (no gorilla or orangutan, ancestor can't be determined)
			if messyToN {
				outputBase = dna.N
			} else {
				outputBase = chimp.Seq[j]
			}
		} else if human.Seq[j] == chimp.Seq[j] { //both descendants match, HCA is same
			outputBase = chimp.Seq[j]
			allMatch++
		} else if (human.Seq[j] == gorilla.Seq[j] || human.Seq[j] == orangutan.Seq[j]) && (human.Seq[j] != dna.N && human.Seq[j] != dna.Gap) { //chimp is different from ancestor, human agrees with gorilla and/or orangutan
			outputBase = human.Seq[j]
			chimpChange++
		} else if (chimp.Seq[j] == gorilla.Seq[j] || human.Seq[j] == gorilla.Seq[j]) && (gorilla.Seq[j] != dna.N && gorilla.Seq[j] != dna.Gap) {
			outputBase = gorilla.Seq[j]
			gorillaVote++
		} else if (human.Seq[j] == orangutan.Seq[j] || chimp.Seq[j] == orangutan.Seq[j]) && (orangutan.Seq[j] != dna.N && orangutan.Seq[j] != dna.Gap) {
			outputBase = orangutan.Seq[j]
			orangutanVote++
		} else {
			if chimp.Seq[j] != dna.Gap && !messyToN {
				outputBase = chimp.Seq[j]
			} else {
				outputBase = dna.N
			}
			messyBase++
		}
		answer.Seq = append(answer.Seq, outputBase)
	}
	return answer
}

//like humanInsertion, determines if the chimp base in the alignment is part of an insertion. Bool return.
func chimpInsertion(records []fasta.Fasta, j int) bool {
	human := records[0]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	if (chimp.Seq[j] != dna.Gap) && (human.Seq[j] == dna.Gap && gorilla.Seq[j] == dna.Gap && orangutan.Seq[j] == dna.Gap) {
		return true
	}
	return false
}
