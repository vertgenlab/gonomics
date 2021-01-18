package gene

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
)

// GtfToGene converts a gtf record into a Gene data structure
// WARNING: If multiple isoforms are present, only the isoform the longest CDS is used
func GtfToGene(g *gtf.Gene, ref []*fasta.Fasta) *Gene {
	answer := new(Gene)
	gtf.MoveCanonicalToZero(g)
	transcript := g.Transcripts[0]

	answer.id = g.GeneID
	answer.symbol = g.GeneName
	answer.posStrand = transcript.Strand

	if transcript.Strand { // If on the positive posStrand
		processExonsPos(transcript, ref, answer)
	} else { // If on the negative posStrand
		processExonsNeg(transcript, ref, answer)
	}

	answer.utrFive.start = 0
	answer.utrThree.end = len(answer.cdnaSeq)
	answer.utrFive.seq = answer.cdnaSeq[answer.utrFive.start:answer.utrFive.end]
	answer.utrThree.seq = answer.cdnaSeq[answer.utrThree.start:answer.utrThree.end]
	answer.codingSeq.start = answer.utrFive.end
	answer.codingSeq.end = answer.utrThree.start
	answer.codingSeq.seq = answer.cdnaSeq[answer.codingSeq.start:answer.codingSeq.end]
	answer.protSeq = dna.TranslateSeq(answer.codingSeq.seq)

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
	answer.orig.utrFive.start = 0
	answer.orig.utrFive.end = answer.utrFive.end
	answer.orig.utrThree.start = answer.utrThree.start
	answer.orig.utrThree.end = answer.utrThree.end
	answer.orig.utrFive.seq = answer.orig.cdnaSeq[answer.utrFive.start:answer.utrFive.end]
	answer.orig.utrThree.seq = answer.orig.cdnaSeq[answer.utrThree.start:answer.utrThree.end]
	answer.orig.codingSeq.start = answer.utrFive.end
	answer.orig.codingSeq.end = answer.utrThree.start
	answer.orig.codingSeq.seq = answer.orig.cdnaSeq[answer.codingSeq.start:answer.codingSeq.end]
	return answer
}

func processExonsPos(transcript *gtf.Transcript, ref []*fasta.Fasta, answer *Gene) {
	// set simple values and initialize slices
	answer.startPos = transcript.Start - 1
	fastaMap := fasta.FastaMap(ref)
	answer.genomeSeq = make([]dna.Base, transcript.End-(transcript.Start-1))
	copy(answer.genomeSeq, fastaMap[transcript.Chr][transcript.Start-1:transcript.End])
	answer.cdnaSeq = make([]dna.Base, 0, len(answer.genomeSeq))
	answer.featureArray = make([]Feature, len(answer.genomeSeq))
	answer.cdsStarts = make([]int, 0, len(transcript.Exons))
	answer.cdsEnds = make([]int, 0, len(transcript.Exons))

	// loop through exons (these are read in order of 5' UTR < CDS < 3' UTR
	var i int
	var prevExonEnd int = answer.startPos
	var currCdsPos int32
	var answerIdxStart, answerIdxEnd int
	for _, val := range transcript.Exons {
		for i = prevExonEnd - answer.startPos; i < val.Start-1-answer.startPos; i++ {
			answer.featureArray[i] = Intron
		}
		prevExonEnd = val.End

		// PROCESS 5' UTR
		if val.FiveUtr != nil {
			answerIdxStart = val.FiveUtr.Start - 1 - answer.startPos // pos relative to answer.genomicSeq
			answerIdxEnd = val.FiveUtr.End - answer.startPos         // pos relative to answer.genomicSeq

			for i = answerIdxStart; i < answerIdxEnd; i++ {
				answer.featureArray[i] = UtrFive
			}

			answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answerIdxStart:answerIdxEnd]...) // get the seq from the ref fasta
			answer.utrFive.end = answerIdxEnd                                                         // answer is in cDNA coordinates (i.e. first base of 5' UTR is zero).
		}

		// PROCESS CDS
		if val.Cds != nil {
			answerIdxStart = val.Cds.Start - 1 - answer.startPos
			answerIdxEnd = val.Cds.End - answer.startPos
			answer.cdsStarts = append(answer.cdsStarts, answerIdxStart)
			answer.cdsEnds = append(answer.cdsEnds, answerIdxEnd-1) // end pos is closed
			answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answerIdxStart:answerIdxEnd]...)

			for i = answerIdxStart; i < answerIdxEnd; i++ {
				answer.featureArray[i] = Feature(currCdsPos)
				currCdsPos++
			}
		}

		// PROCESS 3' UTR
		if val.ThreeUtr != nil {
			answerIdxStart = val.ThreeUtr.Start - 1 - answer.startPos
			answerIdxEnd = val.ThreeUtr.End - answer.startPos

			for i = answerIdxStart; i < answerIdxEnd; i++ {
				answer.featureArray[i] = UtrThree
			}

			// if unset, then define the start
			if answer.utrThree.start == 0 { // by the time the read get to 3' UTR, the start will be defined as the current len of cDNA (5' UTR + CDS)
				answer.utrThree.start = len(answer.cdnaSeq)
			}

			answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answerIdxStart:answerIdxEnd]...)
		}
	}
}

func processExonsNeg(transcript *gtf.Transcript, ref []*fasta.Fasta, answer *Gene) {
	// set simple values and initialize slices
	answer.startPos = transcript.End - 1
	fastaMap := fasta.FastaMap(ref)
	answer.genomeSeq = make([]dna.Base, transcript.End-(transcript.Start-1))
	copy(answer.genomeSeq, fastaMap[transcript.Chr][transcript.Start-1:transcript.End])
	dna.ReverseComplement(answer.genomeSeq)
	answer.cdnaSeq = make([]dna.Base, 0, len(answer.genomeSeq))
	answer.featureArray = make([]Feature, len(answer.genomeSeq))
	answer.cdsStarts = make([]int, 0, len(transcript.Exons))
	answer.cdsEnds = make([]int, 0, len(transcript.Exons))

	// loop through exons backwards (these are read in order of 5' UTR < CDS < 3' UTR
	var i, k int
	var prevExonEnd int = answer.startPos
	var val *gtf.Exon
	var currCdsPos int32
	var answerIdxStart, answerIdxEnd int
	for k = len(transcript.Exons) - 1; k >= 0; k-- {
		val = transcript.Exons[k]

		for i = answer.startPos - prevExonEnd; i < answer.startPos-(val.End-1); i++ {
			answer.featureArray[i] = Intron
		}
		prevExonEnd = val.Start - 2

		// PROCESS 5' UTR
		if val.FiveUtr != nil {
			answerIdxStart = answer.startPos - (val.FiveUtr.End - 1)
			answerIdxEnd = answer.startPos - (val.FiveUtr.Start - 2)

			for i = answerIdxStart; i < answerIdxEnd; i++ {
				answer.featureArray[i] = UtrFive
			}

			answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answerIdxStart:answerIdxEnd]...)
			answer.utrFive.end = answerIdxEnd
		}

		// PROCESS CDS
		if val.Cds != nil {
			answerIdxStart = answer.startPos - (val.Cds.End - 1)
			answerIdxEnd = answer.startPos - (val.Cds.Start - 2)

			answer.cdsStarts = append(answer.cdsStarts, answerIdxStart)
			answer.cdsEnds = append(answer.cdsEnds, answerIdxEnd-1) // end pos is closed
			answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answerIdxStart:answerIdxEnd]...)

			for i = answerIdxStart; i < answerIdxEnd; i++ {
				answer.featureArray[i] = Feature(currCdsPos)
				currCdsPos++
			}
		}

		// PROCESS 3' UTR
		if val.ThreeUtr != nil {
			answerIdxStart = answer.startPos - (val.ThreeUtr.End - 1)
			answerIdxEnd = answer.startPos - (val.ThreeUtr.Start - 2)

			for i = answerIdxStart; i < answerIdxEnd; i++ {
				answer.featureArray[i] = UtrThree
			}

			if answer.utrThree.start == 0 { // if start is unset, then this is the first instance of a 3' UTR
				answer.utrThree.start = len(answer.cdnaSeq) // therefore the utrThreeStart is = curr len of cDNA encoding
			}

			answer.cdnaSeq = append(answer.cdnaSeq, answer.genomeSeq[answerIdxStart:answerIdxEnd]...)
		}
	}
}
