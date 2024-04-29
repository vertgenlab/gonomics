package genomeGraph

import (
	"log"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

func LeftGetTwoBit(n *Node, extension int, position int, seq *dnaTwoBit.TwoBit, ans *dnaTwoBit.TwoBit) *dnaTwoBit.TwoBit {

	basesToTake := numbers.Min(seq.Len+position, extension) - seq.Len

	if basesToTake > 0 {

		frag := dnaTwoBit.GetFrag(n.SeqTwoBit, position-basesToTake, position)
		dnaTwoBit.Cat(ans, frag)
	}
	dnaTwoBit.Cat(ans, seq)
	return ans
}

func TwoBitLocalLeftAlign(alpha *dnaTwoBit.TwoBit, beta *dnaTwoBit.TwoBit, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []cigar.ByteCigar, int, int, int, int) {
	rows, columns := alpha.Len+1, beta.Len+1
	var i, j, routeIdx int

	for i = 0; i < rows; i++ {
		m[i][0] = 0
	}
	for j = 0; j < columns; j++ {
		m[0][j] = 0
	}
	for i = 1; i < rows; i++ {
		for j = 1; j < columns; j++ {
			m[i][j], trace[i][j] = cigar.TraceMatrixExtension(m[i-1][j-1], m[i-1][j-1]+scores[dnaTwoBit.GetBase(alpha, uint(i-1))][dnaTwoBit.GetBase(beta, uint(j-1))], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = alpha.Len, alpha.Len
	var route []cigar.ByteCigar
	//traceback starts in top corner
	for i, j, routeIdx = alpha.Len, beta.Len, 0; m[i][j] > 0; {
		if len(route) == 0 {
			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = cigar.AddCigarByte(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = cigar.AddCigarByte(route, curr)
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
			log.Fatalf("Error: unexpected traceback %c\n", trace[i][j])
		}
		minI = i
		minJ = j
	}
	cigar.ReverseBytesCigar(route)
	return m[alpha.Len][beta.Len], route, minI, alpha.Len, minJ, beta.Len
}

func RightGetTwoBit(n *Node, extension int, position int, seq *dnaTwoBit.TwoBit, ans *dnaTwoBit.TwoBit) *dnaTwoBit.TwoBit {
	// if n == nil || n.SeqTwoBit == nil || seq == nil || ans == nil {
	// 	log.Fatalf("Null pointer provided to RightGetTwoBit")
	// 	return nil
	// }
	basesToTake := numbers.Min(n.SeqTwoBit.Len-position, extension)
	dnaTwoBit.Cat(ans, seq)
	if basesToTake > 0 {
		frag := dnaTwoBit.GetFrag(n.SeqTwoBit, position, position+basesToTake)
		if frag != nil { // Check if the fragment is not nil
			dnaTwoBit.Cat(ans, frag)
		}
	}
	return ans
}

// func TwoBitLocalRightAlign(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []cigar.ByteCigar, int, int, int, int) {
// 	//check if size of alpha is larger than m
// 	var currMax int64
// 	var maxI int
// 	var maxJ int
// 	var i, j, routeIdx int
// 	//setting up the first rows and columns
// 	//seting up the rest of the matrix
// 	for i = 0; i < len(alpha)+1; i++ {
// 		for j = 0; j < len(beta)+1; j++ {
// 			if i == 0 && j == 0 {
// 				m[i][j] = 0
// 			} else if i == 0 {
// 				m[i][j] = m[i][j-1] + gapPen
// 				trace[i][j] = 'I'
// 			} else if j == 0 {
// 				m[i][j] = m[i-1][j] + gapPen
// 				trace[i][j] = 'D'
// 			} else {
// 				m[i][j], trace[i][j] = cigar.TraceMatrixExtension(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
// 			}
// 			if m[i][j] > currMax {
// 				currMax = m[i][j]
// 				maxI = i
// 				maxJ = j
// 			}
// 		}
// 	}
// 	var route []cigar.ByteCigar
// 	//traceback starts in top corner
// 	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
// 		//if route[routeIdx].RunLength == 0 {
// 		if len(route) == 0 {
// 			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
// 			route = append(route, curr)
// 		} else if route[routeIdx].Op == trace[i][j] {
// 			route[routeIdx].RunLen += 1
// 		} else {
// 			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
// 			route = append(route, curr)
// 			routeIdx++
// 		}
// 		switch trace[i][j] {
// 		case '=':
// 			i, j = i-1, j-1
// 		case 'X':
// 			i, j = i-1, j-1
// 		case 'I':
// 			j -= 1
// 		case 'D':
// 			i -= 1
// 		default:
// 			log.Fatalf("Error: unexpected traceback with %c\n", trace[i][j])
// 		}
// 	}
// 	cigar.ReverseBytesCigar(route)
// 	return m[maxI][maxJ], route, 0, maxI, 0, maxJ
// }

// // func getLeftBases(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
// // 	seqLen := len(seq)
// // 	basesToTake := int(math.Min(float64(seqLen+refEnd), float64(extension))) - seqLen
// // 	if basesToTake > 0 {
// // 		ans = append(ans, n.Seq[refEnd-basesToTake:refEnd]...)
// // 	}
// // 	ans = append(ans, seq...)
// // 	return ans
// // }

func scoreTwoBitSeedPart(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func scoreTwoBitSeed(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for ; seed != nil; seed = seed.NextPart {
		score += scoreSeedPart(seed, read, scoreMatrix)
	}
	return score
}
