package genomeGraph

import (
	"log"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func NeedlemanWunsch(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []cigar.Cigar) {
	var i, j, routeIdx int
	for i = 0; i < len(alpha)+1; i++ {
		for j = 0; j < len(beta)+1; j++ {
			if i == 0 && j == 0 {
				trace[i][j] = cigar.Match
			} else if i == 0 {
				m[i][j] = m[i][j-1] + gapPen
				trace[i][j] = cigar.Insertion
			} else if j == 0 {
				m[i][j] = m[i-1][j] + gapPen
				trace[i][j] = cigar.Deletion
			} else {
				m[i][j], trace[i][j] = cigar.TripleMaxTraceExtended(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
		}
	}
	var route []cigar.Cigar
	route = append(route, cigar.Cigar{RunLength: 0, Op: trace[len(alpha)][len(beta)]})
	for i, j, routeIdx = len(alpha)-1, len(beta)-1, 0; i > 0 || j > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			route = append(route, cigar.Cigar{RunLength: 1, Op: trace[i][j]})
			routeIdx++
		}
		switch trace[i][j] {
		case cigar.Equal:
			i, j = i-1, j-1
		case cigar.Mismatch:
			i, j = i-1, j-1
		case cigar.Insertion:
			j -= 1
		case cigar.Deletion:
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
	}
	cigar.ReverseCigar(route)
	return m[len(alpha)-1][len(beta)-1], route
}

// FaSeqToNode is a general function used create a new node based on a target fasta, query fasta and a cigar operation.
// In addition, given two indices, it will update start/end for the subset of bases used to create the new Node.
// TODO: Add logic for correct node name annotation convention.
func FaSeqToNode(target fasta.Fasta, query fasta.Fasta, tStart int, qStart int, cigar align.Cigar, index int) (*Node, int, int) {
	switch cigar.Op {
	case align.ColM:
		curr := &Node{Id: uint32(index), Seq: target.Seq[tStart : tStart+int(cigar.RunLength)]}
		return curr, tStart + int(cigar.RunLength), qStart + int(cigar.RunLength)
	case align.ColI:
		ins := &Node{Id: uint32(index), Seq: query.Seq[qStart : qStart+int(cigar.RunLength)]}
		return ins, tStart, qStart + int(cigar.RunLength)
	case align.ColD:
		del := &Node{Id: uint32(index), Seq: target.Seq[tStart : tStart+int(cigar.RunLength)]}
		return del, tStart + int(cigar.RunLength), qStart
	default:
		log.Fatalf("Error: Did not recognize cigar op %d...\n", cigar.Op)
		return nil, 0, 0
	}
}
