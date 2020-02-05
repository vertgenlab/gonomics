package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
	"time"
)

//TODO: add check to make sure reference names in VCF match fasta reference
func FaToGenomeGraph(ref []*fasta.Fasta, vcfs []*vcf.Vcf) *SimpleGraph {
	gg := NewGraph()
	vcfSplit := vcf.VcfSplit(vcfs, ref)

	if len(vcfSplit) != len(ref) {
		log.Fatal("Slice of vcfs do not equal reference length")
	} else {
		for i := 0; i < len(ref); i++ {
			gg = VcfNodesToGraph(gg, ref[i], vcfSplit[i])
		}
	}
	return gg
}

func VcfNodesToGraph(sg *SimpleGraph, chr *fasta.Fasta, vcfs []*vcf.Vcf) *SimpleGraph {
	var vcfFlag = -1

	var curr *Node
	var prev *Node
	var currMatch *Node = nil
	var lastMatch *Node = nil
	var refAllele *Node
	var altAllele *Node

	var idx int64 = 0
	for i := 0; i < len(vcfs); i++ {
		idx, vcfFlag, curr, prev, currMatch, lastMatch, refAllele, altAllele = NodesToGraph(sg, chr, vcfFlag, vcfs[i], idx, curr, prev, currMatch, lastMatch, refAllele, altAllele)
		if idx > int64(len(chr.Seq)) {
			return sg
		}
	}
	//last node case
	if int64(len(chr.Seq))-idx != 0 {
		lastNode := &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:]}
		AddNode(sg, lastNode)
		if vcfFlag == 1 {
			AddEdge(refAllele, lastNode, 1)
			AddEdge(altAllele, lastNode, 1)
		}
		if vcfFlag == 2 || vcfFlag == 3 {
			AddEdge(lastMatch, lastNode, 0.5)
			AddEdge(prev, lastNode, 1)
		}
	}

	return sg
}

func goGraphSmithWatermanMap(gg *SimpleGraph, read *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: ""}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, seedScore, bestScore int64
	var leftPath, rightPath, bestPath []uint32

	var currScore, maxScore int64 = 0, 0
	for bases := 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}
	ext := int(maxScore/600) + len(read.Seq)

	//seedStart := time.Now()
	var currRead *fastq.Fastq = nil
	var seeds []*SeedDev = findSeedsInMap(seedHash, read, seedLen, stepSize, true)
	revCompRead := fastq.Copy(read)
	fastq.ReverseComplement(revCompRead)
	var revCompSeeds []*SeedDev = findSeedsInMap(seedHash, revCompRead, seedLen, stepSize, false)
	seeds = append(seeds, revCompSeeds...)
	SortSeedDevByLen(seeds)
	//seedEnd := time.Now()
	//seedDuration := seedEnd.Sub(seedStart)
	//log.Printf("Have %d seeds, where best is of length %d, and it took %s\n", len(seeds), seeds[0].Length, seedDuration)

	//swStart := time.Now()
	for i = 0; i < len(seeds) && seedCouldBeBetter(seeds[i], bestScore, maxScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		if seeds[i].PosStrand {
			currRead = read
		} else {
			currRead = revCompRead
		}
		seedScore = scoreSeed(seeds[i], currRead)
		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, ext, currRead.Seq[:seeds[i].QueryStart], m, trace)
		rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart+seeds[i].Length), []uint32{}, ext, currRead.Seq[seeds[i].QueryStart+seeds[i].Length:], m, trace)
		currScore = leftScore + seedScore + rightScore
		if currScore > bestScore {
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 1
			}
			bestPath = CatPaths(AddPath(seeds[i].TargetId, leftPath), rightPath)

			currBest.RName = gg.Nodes[bestPath[0]].Name
			currBest.Pos = int64(minTarget)
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(seeds[i].Length), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(int64(minQuery), int64(len(currRead.Seq)), currBest.Cigar)
			currBest.Extra = PathToString(bestPath, gg)
		}
	}
	//swEnd := time.Now()
	//duration := swEnd.Sub(swStart)
	//log.Printf("sw took %s\n", duration)
	return &currBest
}

func goGraphSmithWaterman(gg *SimpleGraph, read *fastq.Fastq, seedHash [][]*SeedBed, seedLen int, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: ""}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, seedScore, bestScore int64
	var leftPath, rightPath, bestPath []uint32

	var currScore, maxScore int64 = 0, 0
	for bases := 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}
	ext := int(maxScore/600) + len(read.Seq)

	var currRead *fastq.Fastq = nil
	var seeds []*SeedDev = findSeedsInSlice(seedHash, read, seedLen, true)
	revCompRead := fastq.Copy(read)
	fastq.ReverseComplement(revCompRead)
	var revCompSeeds []*SeedDev = findSeedsInSlice(seedHash, revCompRead, seedLen, false)
	seeds = append(seeds, revCompSeeds...)
	SortSeedDevByLen(seeds)

	for i = 0; i < len(seeds) && seedCouldBeBetter(seeds[i], bestScore, maxScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		if seeds[i].PosStrand {
			currRead = read
		} else {
			currRead = revCompRead
		}
		seedScore = scoreSeed(seeds[i], currRead)
		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, ext, currRead.Seq[:seeds[i].QueryStart], m, trace)
		rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart+seeds[i].Length), []uint32{}, ext, currRead.Seq[seeds[i].QueryStart+seeds[i].Length:], m, trace)
		currScore = leftScore + seedScore + rightScore
		if currScore > bestScore {
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 1
			}
			bestPath = CatPaths(AddPath(seeds[i].TargetId, leftPath), rightPath)

			currBest.RName = gg.Nodes[bestPath[0]].Name
			currBest.Pos = int64(minTarget)
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(seeds[i].Length), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(int64(minQuery), int64(len(currRead.Seq)), currBest.Cigar)
			currBest.Extra = PathToString(bestPath, gg)
		}
	}
	return &currBest
}

func NodesToGraph(sg *SimpleGraph, chr *fasta.Fasta, vcfFlag int, v *vcf.Vcf, idx int64, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (int64, int, *Node, *Node, *Node, *Node, *Node, *Node) {
	if chr.Name != v.Chr {
		log.Fatalf("Fasta %s does not match vcf name %s", chr.Name, v.Chr)
	}
	if vcfFlag == -1 && v.Pos == 1 {
		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {
			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)
			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt)}

			AddNode(sg, altAllele)
			vcfFlag = 1
			idx++
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, lastMatch)
			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}
			AddNode(sg, prev)
			AddEdge(lastMatch, prev, 0.5)
			idx = v.Pos
			vcfFlag = 2
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {
			lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Alt)}
			AddNode(sg, lastMatch)
			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Alt[1:len(v.Ref)])}
			AddNode(sg, prev)
			AddEdge(lastMatch, prev, 0.5)

			idx = v.Pos + int64(len(prev.Seq)) - 1
			vcfFlag = 3
		}
	} else {
		if v.Pos-1-idx > 0 {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: chr.Seq[idx:v.Pos]}
			AddNode(sg, currMatch)

			if lastMatch != nil && vcfFlag != 1 {
				AddEdge(lastMatch, currMatch, 0.5)
			}
			if vcfFlag == 1 {
				AddEdge(refAllele, currMatch, 1)
				AddEdge(altAllele, currMatch, 1)
			}
			if vcfFlag == 2 || vcfFlag == 3 {
				AddEdge(prev, currMatch, 1)
			} //AddEdge(lastMatch, currMatch, 0.5)
			vcfFlag = 0
		}
		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {

			if vcfFlag == 0 {
				currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(currMatch, refAllele, 0.5)
				AddEdge(currMatch, altAllele, 0.5)
			} else if vcfFlag == 1 {
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, curr)
				AddEdge(refAllele, curr, 1)
				refAllele = curr

				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, curr)
				AddEdge(altAllele, curr, 1)
				altAllele = curr
			} else if vcfFlag == 2 {

				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)

				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(lastMatch, refAllele, 0.5)
				AddEdge(prev, altAllele, 1)
				curr = refAllele

			} else if vcfFlag == 3 {
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(prev, refAllele, 1)
				AddEdge(lastMatch, altAllele, 0.5)
				curr = refAllele
			} else {
				log.Fatal("Flag was not set up correctly")
			}
			vcfFlag = 1
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}

			if vcfFlag == 0 {
				AddNode(sg, curr)
				AddEdge(currMatch, curr, 0.5)

				vcfFlag = 2
			} else if vcfFlag == 1 {
				AddNode(sg, curr)
				AddEdge(altAllele, curr, 0.5)
				altAllele = curr
				//flag as SNP because we still need to connect two nodes to the next match
				vcfFlag = 1
			} else if vcfFlag == 2 {

				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}
				AddNode(sg, curr)
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 1)

				if int64(len(chr.Seq))-v.Pos == 0 {
					AddEdge(currMatch, curr, 1)
				} else {
					AddEdge(currMatch, curr, 0.5)
				}
				vcfFlag = 2
			} else if vcfFlag == 3 {
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_alt", Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}
				AddNode(sg, curr)
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 1)

				if int64(len(chr.Seq))-v.Pos == 0 {
					AddEdge(currMatch, curr, 1)
				} else {
					AddEdge(currMatch, curr, 0.5)
				}

				vcfFlag = 2
			} else {
				log.Fatal("Flag was not set up correctly")
			}
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {

			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}

			if vcfFlag == 0 {
				AddNode(sg, curr)
				AddEdge(currMatch, curr, 0.5)
				vcfFlag = 3
			} else if vcfFlag == 1 {
				AddNode(sg, curr)
				AddEdge(refAllele, curr, 0.5)
				refAllele = curr
				vcfFlag = 1
			} else if vcfFlag == 2 {

				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}
				AddNode(sg, curr)
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 1)

				if int64(len(chr.Seq))-v.Pos-1 == 0 {
					AddEdge(currMatch, curr, 1)
				} else {
					AddEdge(currMatch, curr, 0.5)
				}
				vcfFlag = 3
			} else if vcfFlag == 3 {
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name + "_ref", Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}
				AddNode(sg, curr)
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 1)

				if int64(len(chr.Seq))-v.Pos == 0 {
					AddEdge(currMatch, curr, 1)
				} else {
					AddEdge(currMatch, curr, 0.5)
				}
				vcfFlag = 3
			} else {
				log.Fatal("Flag was not set up correctly")
			}

			idx = v.Pos + int64(len(curr.Seq))
		}
		prev = curr
		lastMatch = currMatch

	}

	return idx, vcfFlag, curr, prev, currMatch, lastMatch, refAllele, altAllele

}

func wrapMap(ref *SimpleGraph, r *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, c chan *sam.SamAln) {
	var mappedRead *sam.SamAln
	m, trace := swMatrixSetup(10000)
	mappedRead = goGraphSmithWatermanMap(ref, r, seedHash, seedLen, stepSize, m, trace)
	c <- mappedRead
	//log.Printf("%s\n", sam.SamAlnToString(mappedRead))
}

func wrap(ref *SimpleGraph, r *fastq.Fastq, seedHash [][]*SeedBed, seedLen int, c chan *sam.SamAln) {
	var mappedRead *sam.SamAln
	m, trace := swMatrixSetup(10000)
	mappedRead = goGraphSmithWaterman(ref, r, seedHash, seedLen, m, trace)
	c <- mappedRead
	//log.Printf("%s\n", sam.SamAlnToString(mappedRead))
}

func wrapNoChanMap(ref *SimpleGraph, r *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int) {
	var mappedRead *sam.SamAln
	startOne := time.Now()
	m, trace := swMatrixSetup(10000)
	stopOne := time.Now()
	durationOne := stopOne.Sub(startOne)
	startTwo := time.Now()
	mappedRead = goGraphSmithWatermanMap(ref, r, seedHash, seedLen, stepSize, m, trace)
	stopTwo := time.Now()
	durationTwo := stopTwo.Sub(startTwo)
	log.Printf("%s %s %s\n", sam.SamAlnToString(mappedRead), durationOne, durationTwo)
}

func wrapNoChan(ref *SimpleGraph, r *fastq.Fastq, seedHash [][]*SeedBed, seedLen int) {
	var mappedRead *sam.SamAln
	m, trace := swMatrixSetup(10000)
	mappedRead = goGraphSmithWaterman(ref, r, seedHash, seedLen, m, trace)
	log.Printf("%s\n", sam.SamAlnToString(mappedRead))
}

