package gene

import "github.com/vertgenlab/gonomics/dna"

// Reset reverts all mutations done to a Gene.
func Reset(g *Gene) {
	var hasIndel bool
	var err error
	for _, val := range g.changeLog {
		if len(val.added) != 1 || len(val.removed) != 1 {
			hasIndel = true
			break
		}
	}

	if !hasIndel { // Point mutations can be quickly reversed by performing the inverse function
		for i := len(g.changeLog) - 1; i >= 0; i-- {
			if len(g.changeLog[i].added) == 1 && len(g.changeLog[i].removed) == 1 {
				_, err = Substitution(g, g.changeLog[i].genomePos, g.changeLog[i].removed[0])
				g.changeLog = g.changeLog[:len(g.changeLog)-2]
				if err != nil {
					hasIndel = true
					break
				}
			}
		}
	}

	if hasIndel { // Indels can remove important feature context which is non-trivial to store in the changelog, so indels must restore from backup
		g.startPos = g.orig.startPos
		g.cdsStarts = g.cdsStarts[:len(g.orig.cdsStarts)]
		copy(g.cdsStarts, g.orig.cdsStarts)
		g.cdsEnds = g.cdsEnds[:len(g.orig.cdsEnds)]
		copy(g.cdsEnds, g.orig.cdsEnds)
		g.genomeSeq = g.genomeSeq[:len(g.orig.genomeSeq)]
		copy(g.genomeSeq, g.orig.genomeSeq)
		g.cdnaSeq = g.cdnaSeq[:len(g.orig.cdnaSeq)]
		copy(g.cdnaSeq, g.orig.cdnaSeq)
		g.featureArray = g.featureArray[:len(g.orig.featureArray)]
		copy(g.featureArray, g.orig.featureArray)
		g.codingSeq.start = g.orig.codingSeq.start
		g.codingSeq.end = g.orig.codingSeq.end
		g.codingSeq.seq = g.cdnaSeq[g.codingSeq.start:g.codingSeq.end]
		g.utrFive.start = g.orig.utrFive.start
		g.utrFive.end = g.orig.utrFive.end
		g.utrFive.seq = g.cdnaSeq[g.utrFive.start:g.utrFive.end]
		g.utrThree.start = g.orig.utrThree.start
		g.utrThree.end = g.orig.utrThree.end
		g.utrThree.seq = g.cdnaSeq[g.utrThree.start:g.utrThree.end]
	}

	g.changeLog = nil
	g.protSeq = dna.TranslateSeq(g.codingSeq.seq)
}
