package gene

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
)

// GtfToGoGene converts a gtf record into a Gene data structure
// WARNING: If multiple isoforms are present, only the isoform the longest CDS is used
func GtfToGoGene(g *gtf.Gene, ref []*fasta.Fasta) *Gene {
	answer := new(Gene)
	gtf.MoveCanonicalToZero(g)
	transcript := g.Transcripts[0]

	answer.id = g.GeneID
	answer.posStrand = transcript.Strand

	if transcript.Strand { // If on the positive posStrand
		answer.startPos = transcript.Start - 1

		fastaMap := fasta.FastaMap(ref)
		answer.genomeSeq = make([]dna.Base, transcript.End-(transcript.Start-1))
		copy(answer.genomeSeq, fastaMap[transcript.Chr][transcript.Start-1:transcript.End])
		answer.cdnaSeq = make([]dna.Base, 0, gtf.CdsLength(transcript))
		answer.featureArray = make([]Feature, len(answer.genomeSeq))
		answer.cdsStarts = make([]int, 0, len(transcript.Exons))
		answer.cdsEnds = make([]int, 0, len(transcript.Exons))

		var i int
		var prevExonEnd int = answer.startPos
		var currCdsPos int32
		for _, val := range transcript.Exons {
			for i = prevExonEnd - answer.startPos; i < val.Start-1-answer.startPos; i++ {
				answer.featureArray[i] = Intron
			}
			prevExonEnd = val.End

			if val.FiveUtr != nil {
				for i = val.FiveUtr.Start - 1 - answer.startPos; i < val.FiveUtr.End-answer.startPos; i++ {
					answer.featureArray[i] = UtrFive
				}
			}

			if val.ThreeUtr != nil {
				for i = val.ThreeUtr.Start - 1 - answer.startPos; i < val.ThreeUtr.End-answer.startPos; i++ {
					answer.featureArray[i] = UtrThree
				}
			}

			if val.Cds != nil {
				answer.cdsStarts = append(answer.cdsStarts, val.Cds.Start-1-answer.startPos)
				answer.cdsEnds = append(answer.cdsEnds, val.Cds.End-1-answer.startPos)
				answer.cdnaSeq = append(answer.cdnaSeq, fastaMap[transcript.Chr][val.Cds.Start-1:val.Cds.End]...)

				for i = val.Cds.Start - 1 - answer.startPos; i < val.Cds.End-answer.startPos; i++ {
					answer.featureArray[i] = Feature(currCdsPos)
					currCdsPos++
				}
			}
		}

	} else { // If on the negative posStrand
		answer.startPos = transcript.End - 1

		fastaMap := fasta.FastaMap(ref)
		answer.genomeSeq = make([]dna.Base, transcript.End-(transcript.Start-1))
		copy(answer.genomeSeq, fastaMap[transcript.Chr][transcript.Start-1:transcript.End])
		dna.ReverseComplement(answer.genomeSeq)
		answer.cdnaSeq = make([]dna.Base, 0, gtf.CdsLength(transcript))
		answer.featureArray = make([]Feature, len(answer.genomeSeq))
		answer.cdsStarts = make([]int, 0, len(transcript.Exons))
		answer.cdsEnds = make([]int, 0, len(transcript.Exons))

		var i, k int
		var prevExonEnd int = answer.startPos
		var val *gtf.Exon
		var currCdsPos int32
		for k = len(transcript.Exons) - 1; k >= 0; k-- {
			val = transcript.Exons[k]
			for i = answer.startPos - prevExonEnd; i < answer.startPos-(val.End-1); i++ {
				answer.featureArray[i] = Intron
			}
			prevExonEnd = val.Start - 2

			if val.FiveUtr != nil {
				for i = answer.startPos - (val.FiveUtr.End - 1); i < answer.startPos-(val.FiveUtr.Start-2); i++ {
					answer.featureArray[i] = UtrFive
				}
			}

			if val.ThreeUtr != nil {
				for i = answer.startPos - (val.ThreeUtr.End - 1); i < answer.startPos-(val.ThreeUtr.Start-2); i++ {
					answer.featureArray[i] = UtrThree
				}
			}

			if val.Cds != nil {
				answer.cdsStarts = append(answer.cdsStarts, answer.startPos-(val.Cds.End-1))
				answer.cdsEnds = append(answer.cdsEnds, answer.startPos-(val.Cds.Start-1))
				answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answer.startPos-(val.Cds.End-1):answer.startPos-(val.Cds.Start-2)]...)

				for i = answer.startPos - (val.Cds.End - 1); i < answer.startPos-(val.Cds.Start-2); i++ {
					answer.featureArray[i] = Feature(currCdsPos)
					currCdsPos++
				}
			}
		}
	}
	answer.orig.startPos = answer.startPos
	answer.orig.cdsStarts = make([]int, len(answer.cdsStarts))
	copy(answer.orig.cdsStarts, answer.cdsStarts)
	answer.orig.cdsEnds = make([]int, len(answer.cdsEnds))
	copy(answer.orig.cdsEnds, answer.cdsEnds)
	answer.orig.genomeSeq = make([]dna.Base, len(answer.genomeSeq))
	copy(answer.orig.genomeSeq, answer.genomeSeq)
	answer.orig.cdnaSeq = make([]dna.Base, len(answer.cdnaSeq))
	copy(answer.orig.cdnaSeq, answer.cdnaSeq)
	answer.orig.featureArray = make([]Feature, len(answer.featureArray))
	copy(answer.orig.featureArray, answer.featureArray)
	return answer
}
