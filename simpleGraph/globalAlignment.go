package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func nmMatrixSetup(size int64) ([][]int64, [][]rune) {
	m := make([][]int64, size)
	trace := make([][]rune, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]rune, size)
	}
	return m, trace
}

func NeedlemanWunsch(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]rune) (int64, []*cigar.Cigar) {
	var i, j, routeIdx int
	for i = 0; i < len(alpha)+1; i++ {
		for j = 0; j < len(beta)+1; j++ {
			if i == 0 && j == 0 {
				m[i][j] = 'M'
			} else if i == 0 {
				m[i][j] = m[i][j-1] + gapPen
				trace[i][j] = 'I'
			} else if j == 0 {
				m[i][j] = m[i-1][j] + gapPen
				trace[i][j] = 'D'
			} else {
				m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
		}
	}
	var route []*cigar.Cigar
	route = append(route, &cigar.Cigar{RunLength: 0, Op: trace[len(alpha)][len(beta)]})
	for i, j, routeIdx = len(alpha)-1, len(beta)-1, 0; i > 0 || j > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			route = append(route, &cigar.Cigar{RunLength: 1, Op: trace[i][j]})
			routeIdx++
		}
		switch trace[i][j] {
		case '=':
			i, j = i-1, j-1
		case 'X':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
	}
	reverseCigarPointer(route)
	return m[len(alpha)-1][len(beta)-1], route
}

// FaSeqToNode is a general function used create a new node based on a target fasta, query fasta and a cigar operation.
// In addition, given two indices, it will update start/end for the subset of bases used to create the new Node.
// TODO: Add logic for correct node name annotation convention
func FaSeqToNode(target *fasta.Fasta, query *fasta.Fasta, tStart int, qStart int, cigar align.Cigar, index int) (*Node, int, int) {
	switch cigar.Op {
	case align.ColM:
		curr := Node{Id: uint32(index), Name: fmt.Sprintf("%s_%d_%d_%s_%d_%d", target.Name, tStart, tStart+int(cigar.RunLength), query.Name, qStart, qStart+int(cigar.RunLength)), Seq: target.Seq[tStart : tStart+int(cigar.RunLength)]}
		return &curr, tStart + int(cigar.RunLength), qStart + int(cigar.RunLength)
	case align.ColI:
		ins := Node{Id: uint32(index), Name: fmt.Sprintf("%s_%d_%d_%s_%d_%d", target.Name, tStart, tStart, query.Name, qStart, qStart+int(cigar.RunLength)), Seq: query.Seq[qStart : qStart+int(cigar.RunLength)]}
		return &ins, tStart, qStart + int(cigar.RunLength)
	case align.ColD:
		del := Node{Id: uint32(index), Name: fmt.Sprintf("%s_%d_%d_%s_%d_%d", target.Name, tStart, tStart+int(cigar.RunLength), query.Name, qStart, qStart), Seq: target.Seq[tStart : tStart+int(cigar.RunLength)]}
		return &del, tStart + int(cigar.RunLength), qStart
	default:
		log.Fatalf("Error: Did not recognize cigar op %d...\n", cigar.Op)
	}
	return nil, 0, 0
}
