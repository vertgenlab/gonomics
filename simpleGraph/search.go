package simpleGraph

import (
	"fmt"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//"sync"
)

func VariantGraph(reference []*fasta.Fasta, vcfs []*vcf.Vcf) *SimpleGraph {
	gg := NewGraph()
	vcfSplit := vcf.VcfSplit(vcfs, reference)
	if len(vcfSplit) != len(reference) {
		log.Fatal("Slice of vcfs do not equal reference length")
	} else {
		for i := 0; i < len(reference); i++ {
			gg = vChrGraph(gg, reference[i], vcfSplit[i])
		}
	}
	return gg
}

func vChrGraph(genome *SimpleGraph, chr *fasta.Fasta, vcfsChr []*vcf.Vcf) *SimpleGraph {
	fasta.ToUpper(chr)
	if len(vcfsChr) == 0 {
		log.Fatalf("Error: vcf file is empty, please try again...\n")
	}
	var currMatch *Node = nil
	var lastMatch *Node = nil
	var refAllele, altAllele *Node
	var prev []*Node = nil
	var weight float32 = 0
	var i, edge int //, j, edge int
	var j int
	var index int64 = 0
	//for debuging, max index can go up to in the for loop:
	//var lastPos int64 = int64(len(chr.Seq)) - vcfsChr[len(vcfsChr)-1].Pos - 1
	//var lastV *vcf.Vcf = *vcf.Vcf{Pos: 0}
	for i = 0; i < len(vcfsChr); i++ {
		//trivial case
		if vcfsChr[i].Pos-index > 0 {
			currMatch = &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d", chr.Name, index+1), Seq: chr.Seq[index : vcfsChr[i].Pos-1], Prev: nil, Next: make([]*Edge, 0, 2)}
			AddNode(genome, currMatch)
			//	currMatch = &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d", chr.Name, index+1), Seq: []dna.Base{chr.Seq[index]}, Prev: nil, Next: make([]*Edge, 0, 2)}
			//}
		} else {
			log.Printf("Warning: there are two vcf records in %s, postion %d...\n", chr.Name, vcfsChr[i].Pos)
			fmt.Printf("%s\t%s\t%d\t%s\t%s\n%s\t%s\t%d\t%s\t%s\n", vcfsChr[i-1].Format, vcfsChr[i-1].Chr, vcfsChr[i-1].Pos, vcfsChr[i-1].Ref, vcfsChr[i-1].Alt, vcfsChr[i].Format, vcfsChr[i].Chr, vcfsChr[i].Pos, vcfsChr[i].Ref, vcfsChr[i].Alt)
			continue
		}
		if lastMatch != nil {
			AddEdge(lastMatch, currMatch, 0.5)
		}
		//if prev != nil {
		//	weight = float32(1) / float32(len(prev))
		//	for edge = 0; edge < len(prev); edge++ {
		//		AddEdge(prev[edge], currMatch, weight)
		//	}
		//}

		prev = make([]*Node, 0, 2)

		if isSNP(vcfsChr[i]) {
			refAllele = &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d_SNP", chr.Name, vcfsChr[i].Pos), Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil}
			AddNode(genome, refAllele)
			AddEdge(currMatch, refAllele, 0.5)

			altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d_SNP", chr.Name, vcfsChr[i].Pos), Seq: dna.StringToBases(vcfsChr[i].Alt), Prev: nil, Next: nil}
			AddNode(genome, altAllele)
			AddEdge(currMatch, altAllele, 0.5)

			prev = append(prev, refAllele)
			prev = append(prev, altAllele)
			index = vcfsChr[i].Pos
			for j = i + 1; j < len(vcfsChr)-1; j++ {
				if isSNP(vcfsChr[j-1]) && isSNP(vcfsChr[j]) && vcfsChr[j].Pos-1 == vcfsChr[j-1].Pos {
					refAllele.Seq = append(refAllele.Seq, dna.StringToBases(vcfsChr[j].Ref)...)
					altAllele.Seq = append(altAllele.Seq, dna.StringToBases(vcfsChr[j].Alt)...)
					index = vcfsChr[j].Pos
				} else {
					lastMatch = currMatch
					i = j - 1
					break
				}
			}
		}
		if isINS(vcfsChr[i]) {
			currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
			insertion := &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d_INS", chr.Name, vcfsChr[i].Pos), Seq: dna.StringToBases(vcfsChr[i].Alt)[1:], Prev: nil, Next: nil}
			AddNode(genome, insertion)
			AddEdge(currMatch, insertion, 0.5)
			prev = append(prev, insertion)
			index = vcfsChr[i].Pos
		}
		if isDEL(vcfsChr[i]) {
			currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Alt)...)
			deletion := &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d_DEL", chr.Name, index+1), Seq: dna.StringToBases(vcfsChr[i].Ref)[1:], Prev: nil, Next: nil}
			//deletion.Name += fmt.Sprintf("_del_%d_%d", vcfsChr[i].Pos+1, vcfsChr[i].Pos+int64(len(deletion.Seq)))
			AddNode(genome, deletion)
			AddEdge(currMatch, deletion, 0.5)
			prev = append(prev, deletion)
			index = vcfsChr[i].Pos + int64(len(deletion.Seq))
		}
		lastMatch = currMatch
		//lastV = vcfsChr[i]
		//currMatch.Name+= fmt.Sprintf("_%d", len(currMatch.Seq))

	}
	//Case: last node
	lastNode := &Node{Id: uint32(len(genome.Nodes)), Name: fmt.Sprintf("%s_%d_%d", chr.Name, index+1, int64(len(chr.Seq))-index), Seq: chr.Seq[index:], Prev: make([]*Edge, 0, len(prev)), Next: nil}
	AddNode(genome, lastNode)
	weight = float32(1) / float32(len(prev))
	for edge = 0; edge < len(prev); edge++ {
		AddEdge(prev[edge], lastNode, weight)
	}
	return genome
}

func AlignTraversalFwd(rightNode *Node, seq []dna.Base, start int, currentPath []uint32, extention int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, []uint32) {
	currentPath = AddPath(rightNode.Id, currentPath)
	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	if len(seq) >= extention {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extension.\n")
	}
	var availableBases int = len(rightNode.Seq) - start + len(seq)

	var targetLength int = common.Min(availableBases, extention)
	var basesToTake int = targetLength - len(seq)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(rightNode.Seq), start, targetLength, basesToTake)
	var s []dna.Base = make([]dna.Base, targetLength)

	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], rightNode.Seq[start:start+basesToTake])

	/*
		var base int
		for base = 0; base < len(seq); base++ {
			s[base] = seq[base]
		}
		for base = 0; base < len(n.Seq[start:start+basesToTake]); base++ {
			s[len(seq)+base] = n.Seq[start+base]
		}*/
	if availableBases >= extention || len(rightNode.Next) == 0 {
		score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, queryEnd, currentPath
	} else {

		for _, i := range rightNode.Next {
			bestScore = -1
			alignment, score, queryEnd, path = AlignTraversalFwd(i.Dest, s, 0, currentPath, extention, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestQueryEnd, bestPath
	}
}

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, extention int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, int, []uint32) {
	currentPath = AddPath(n.Id, currentPath)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, extention)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)
	/*var base int
	for base = 0; base < targetLength; base++ {
		if base < basesToTake {
			s[base] = n.Seq[refEnd-basesToTake+base]
		} else {
			s[base] = seq[base]
		}
	}*/
	if availableBases >= extention || len(n.Next) == 0 {
		score, alignment, refStart, _, queryStart, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		//tmp := make([]uint32, len(currentPath))
		//copy(tmp, currentPath)
		for _, i := range n.Prev {
			alignment, score, refStart, queryStart, path = AlignReverseGraphTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, extention, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
		reversePath(bestPath)
		return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
	}
}

func PointToBases(currSeq []dna.Base, incomingSeq []dna.Base, currStart int, incomingStart int, currEnd int, incomingEnd int, currFront bool) {
	newSize := currEnd + incomingEnd - currStart - incomingStart
	var i int
	if len(currSeq) < newSize {
		currSeq = append(currSeq, make([]dna.Base, newSize-len(currSeq))...)
	}
	if currFront {
		for i = currEnd; i < len(incomingSeq[incomingStart:incomingEnd]); i++ {
			currSeq[i] = incomingSeq[incomingStart+i]
		}
	}
	if !currFront {
		for i = 0; i < len(currSeq); i++ {
			currSeq[incomingEnd-incomingStart-1] = currSeq[i]
		}
		for i = 0; i < len(incomingSeq[incomingStart:incomingEnd]); i++ {
			currSeq[i] = incomingSeq[incomingStart+i]
		}
	}
}
