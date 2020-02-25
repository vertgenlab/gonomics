package simpleGraph

/*
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
	fasta.ToUpper(chr)
	var vcfFlag = -1
	var curr *Node
	var prev *Node
	var currMatch *Node = nil
	var lastMatch *Node = nil
	var refAllele *Node
	var altAllele *Node

	var idx int64 = 0
	var i int
	for i = 0; i < len(vcfs); i++ {
		idx, vcfFlag, curr, prev, currMatch, lastMatch, refAllele, altAllele = NodesToGraph(sg, chr, vcfFlag, vcfs[i], idx, curr, prev, currMatch, lastMatch, refAllele, altAllele)
		if i < len(vcfs) -1 {
			if isSNP(vcfs[i+1]) && isSNP(vcfs[i]) && vcfs[i].Pos+1 == vcfs[i+1].Pos {
				refAllele.Seq = append(refAllele.Seq, dna.StringToBases(vcfs[i+1].Ref)...)
				altAllele.Seq = append(altAllele.Seq, dna.StringToBases(vcfs[i+1].Alt)...)
				i++
				idx++
			}
		}

		if idx > int64(len(chr.Seq)) {
			return sg
		}
	}
	//last node case
	if int64(len(chr.Seq))-idx != 0 {
		lastNode := &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:]}
		AddNode(sg, lastNode)
		if vcfFlag == 1 {
			lastNode.Seq = lastNode.Seq[1:]
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

func NodeSplitByNs(sg *SimpleGraph, currMatch *Node, chr *fasta.Fasta, index int64, end int64) *Node {
	var inRegion bool = false
	var start int64 = 0
	for ; index < end; index++ {
		if dna.DefineBase(chr.Seq[start]) && inRegion == false {
			inRegion = false
			currMatch.Seq = chr.Seq[start:index]
			start = index
		} else if dna.DefineBase(chr.Seq[start]) && inRegion == true {
			newMatch := &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			AddNode(sg, newMatch)
			AddEdge(currMatch, newMatch, 1)
			inRegion = true
			currMatch = newMatch
		}

	}

	//currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
	return currMatch
}

type noNsBed struct {
	Start int32
	End   int32
}

func findUngap(fa *fasta.Fasta, index int64, end int64) []*noNsBed {
	var answer []*noNsBed
	var inRegion bool = false
	var startIndex int32 = 0
	for ; index < end; index++ {
		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
			inRegion = true
			startIndex = int32(index)
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
			answer = append(answer, &noNsBed{Start: startIndex, End: int32(index)})
			inRegion = false
		}

	}
	if inRegion == true {
		answer = append(answer, &noNsBed{Start: int32(startIndex), End: int32(index)})
	}
	return answer
}

func NodesToGraph(sg *SimpleGraph, chr *fasta.Fasta, vcfFlag int, v *vcf.Vcf, idx int64, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (int64, int, *Node, *Node, *Node, *Node, *Node, *Node) {
	if chr.Name != v.Chr {
		log.Fatalf("Fasta %s does not match vcf name %s", chr.Name, v.Chr)
	}
	/*currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:v.Pos-1]}
	AddNode(sg, currMatch)
	if lastMatch != nil && vcfFlag != 1 {
		AddEdge(lastMatch, currMatch, 0.5)
	}
	if vcfFlag == 2 || vcfFlag == 3 {
		AddEdge(prev, currMatch, 1)
	}
	vcfFlag = 0
	//log.Printf("Nodes that dont have sequence: %v, =%d, vcf.pos=%d\n", chr.Seq[idx:v.Pos], idx, v.Pos)

	noN := findUngap(chr, idx, v.Pos)
	if len(noN) == 0 || idx >= v.Pos {
		currMatch = lastMatch
	} else {
		var newMatch *Node = &Node{}
		for i := 0; i < len(noN); i++ {
			if i == 0 {
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[noN[i].Start:noN[i].End]}
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
				}
			} else {
				//AddEdge(currMatch, gappedNode, 1)
				//gappedNode := &Node{Name: "Gap", Seq: dna.StringToBases("NNNNNNNNNNNNNNNNNNNNNN")}
				//AddEdge(currMatch, gappedNode, 1)
				newMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[noN[i].Start:noN[i].End]}
				AddNode(sg, newMatch)
				//AddEdge(currMatch, newMatch, 1)
				currMatch = newMatch
			}
		}
		vcfFlag = 0
	}

	if isSNP(v) {
		//if currMatch != nil || {
		//	currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
		//}
		refAllele, altAllele = SnpToNodesEdges(sg, v, vcfFlag, chr, curr, prev, currMatch, lastMatch, refAllele, altAllele)
		vcfFlag = 1
		idx = v.Pos
	}
	if isINS(v) {
		curr, vcfFlag = InsToNodesEdges(sg, v, vcfFlag, chr, curr, prev, currMatch, lastMatch, refAllele, altAllele)
		idx = v.Pos
	}
	if isDEL(v) {
		curr, vcfFlag = DelToNodesEdges(sg, v, vcfFlag, chr, curr, prev, currMatch, lastMatch, refAllele, altAllele)
		idx = v.Pos + int64(len(curr.Seq))
	}
	prev = curr
	lastMatch = currMatch
	return idx, vcfFlag, curr, prev, currMatch, lastMatch, refAllele, altAllele
}

func NodesToGraph(sg *SimpleGraph, chr *fasta.Fasta, vcfFlag int, v *vcf.Vcf, idx int64, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (int64, int, *Node, *Node, *Node, *Node, *Node, *Node) {
	if chr.Name != v.Chr {
		log.Fatalf("Fasta %s does not match vcf name %s", chr.Name, v.Chr)
	}

	if vcfFlag == -1 {

		currMatch = NodeSplitByNs(sg, lastMatch, chr, idx, v.Pos)

		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {
			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)
			AddEdge(currMatch, refAllele, 0.5)
			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddEdge(currMatch, altAllele, 0.5)
			AddNode(sg, altAllele)
			vcfFlag = 1
			idx = v.Pos - 1
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			//lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			//AddNode(sg, lastMatch)
			//lastMatch = NodeSplitByNs(sg, lastMatch, chr.Seq, idx, v.Pos)

			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			prev.Seq = prev.Seq[1:]
			AddNode(sg, prev)
			AddEdge(currMatch, prev, 0.5)
			idx = v.Pos
			vcfFlag = 2
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {
			//lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			//AddNode(sg, lastMatch)
			//lastMatch = NodeSplitByNs(sg, lastMatch, chr.Seq, idx, v.Pos)
			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, prev)
			AddEdge(currMatch, prev, 0.5)

			idx = v.Pos + int64(len(prev.Seq)) - 1
			lastMatch = currMatch
			vcfFlag = 3
		}
	} else {
		if v.Pos-1-idx > 0 {
			//currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			if lastMatch != nil && vcfFlag != 1 {
				AddEdge(lastMatch, currMatch, 0.5)
			}
			if vcfFlag == 1 {
				AddEdge(refAllele, currMatch, 1)
				AddEdge(altAllele, currMatch, 1)
			}
			if vcfFlag == 2 || vcfFlag == 3 {
				AddEdge(prev, currMatch, 1)
			}

			currMatch = NodeSplitByNs(sg, lastMatch, chr, idx, v.Pos)

			vcfFlag = 0
		}
		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {

			if vcfFlag == 0 {
				currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				AddEdge(currMatch, refAllele, 0.5)

				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)
				AddEdge(currMatch, altAllele, 0.5)

			} else if vcfFlag == 1 {
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, curr)
				AddEdge(refAllele, curr, 1)
				refAllele = curr

				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, curr)
				AddEdge(altAllele, curr, 1)
				altAllele = curr
			} else if vcfFlag == 2 {

				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)

				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(lastMatch, refAllele, 0.5)
				AddEdge(prev, altAllele, 1)
				curr = refAllele

			} else if vcfFlag == 3 {
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(prev, refAllele, 1)
				AddEdge(lastMatch, altAllele, 0.5)
				curr = refAllele
			} else {
				log.Fatal("Flag was not set up correctly")
			}
			lastMatch = currMatch
			vcfFlag = 1
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			curr.Seq= curr.Seq[1:]
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

				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
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
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
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
			prev = curr
			lastMatch = currMatch
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {

			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)[1:]}

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

				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
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
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
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

			idx = v.Pos + int64(len(curr.Seq)) - 1
		}
		prev = curr
		lastMatch = currMatch

	}

	return idx, vcfFlag, curr, prev, currMatch, lastMatch, refAllele, altAllele

}

func SnpToNodesEdges(sg *SimpleGraph, snp *vcf.Vcf, lastNodeFlag int, chr *fasta.Fasta, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (*Node, *Node) {
	//Special case if last variant called was a SNP
	if lastNodeFlag == 1 {
		refAllele.Seq = append(refAllele.Seq, dna.StringToBases(snp.Ref)...)
		altAllele.Seq = append(altAllele.Seq, dna.StringToBases(snp.Alt)...)
		//curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Ref)}
		//AddNode(sg, curr)
		//AddEdge(refAllele, curr, 1)
		//refAllele = curr

		//curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Alt)}
		//AddNode(sg, curr)
		//AddEdge(altAllele, curr, 1)
		//altAllele = curr
	} else {
		refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Ref)}
		AddNode(sg, refAllele)
		altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Alt)}
		AddNode(sg, altAllele)
		//last node was a match
		if lastNodeFlag == 0 {
			if currMatch != nil {
				currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
			}

			AddEdge(currMatch, refAllele, 0.5)
			AddEdge(currMatch, altAllele, 0.5)
		}//last node was an insertion
		if lastNodeFlag == 2 {
			AddEdge(lastMatch, refAllele, 0.5)
			AddEdge(prev, altAllele, 1)
		}//last node was a deletion
		if lastNodeFlag == 3 {
			AddEdge(prev, refAllele, 1)
			AddEdge(lastMatch, altAllele, 0.5)
		}
	}
	return refAllele, altAllele
}

func InsToNodesEdges(sg *SimpleGraph, ins *vcf.Vcf, vcfFlag int, chr *fasta.Fasta, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (*Node, int) {
	if vcfFlag == 0 || vcfFlag == 1 {
		curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(ins.Alt[1:])}
		AddNode(sg, curr)
		if vcfFlag == 0 {
			AddEdge(currMatch, curr, 0.5)
			vcfFlag = 2
		} else {
			//flag as SNP because we still need to connect two nodes to the next match
			AddEdge(altAllele, curr, 0.5)
			altAllele = curr
			vcfFlag = 1
		}
	} else if vcfFlag == 2 || vcfFlag == 3 {
		if vcfFlag == 2 {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(ins.Ref)}
			curr = &Node{Id: uint32(len(sg.Nodes)+1), Name: chr.Name, Seq: dna.StringToBases(ins.Alt[1:])}
		} else {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(ins.Alt)}
			curr = &Node{Id: uint32(len(sg.Nodes)+1), Name: chr.Name, Seq: dna.StringToBases(ins.Ref[1:])}
		}
		AddNode(sg, currMatch)
		AddNode(sg, curr)
		AddEdge(lastMatch, currMatch, 0.5)
		AddEdge(prev, currMatch, 1)
		if int64(len(chr.Seq))-ins.Pos == 0 {
			AddEdge(currMatch, curr, 1)
		} else {
			AddEdge(currMatch, curr, 0.5)
		}
		vcfFlag = 2
	}
	return curr, vcfFlag
}

func DelToNodesEdges(sg *SimpleGraph, del *vcf.Vcf, vcfFlag int, chr *fasta.Fasta, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (*Node, int) {
	if vcfFlag == 0 || vcfFlag == 1 {
		curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(del.Ref[1:])}
		AddNode(sg, curr)
		if vcfFlag == 0 {
			AddEdge(currMatch, curr, 0.5)
			vcfFlag = 3
		} else {
			//flag as SNP because we still need to connect two nodes to the next match
			AddEdge(refAllele, curr, 0.5)
			altAllele = curr
			vcfFlag = 1
		}
	} else if vcfFlag == 2 || vcfFlag == 3 {
		if vcfFlag == 2 {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(del.Alt)}
			curr = &Node{Id: uint32(len(sg.Nodes)+1), Name: chr.Name, Seq: dna.StringToBases(del.Ref[1:])}
		} else {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(del.Ref)}
			curr = &Node{Id: uint32(len(sg.Nodes)+1), Name: chr.Name, Seq: dna.StringToBases(del.Alt[1:])}
		}
		AddNode(sg, currMatch)
		AddNode(sg, curr)
		AddEdge(lastMatch, currMatch, 0.5)
		AddEdge(prev, currMatch, 1)
		if int64(len(chr.Seq))-del.Pos == 0 {
			AddEdge(currMatch, curr, 1)
		} else {
			AddEdge(currMatch, curr, 0.5)
		}
		vcfFlag = 3
	}
	return curr, vcfFlag
}

func isSNP(v *vcf.Vcf) bool {
	var truth bool = false
	if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {
		truth = true
	}
	return truth
}

func createINS(sg *SimpleGraph, v *vcf.Vcf, chr string) *Node {
	curr := Node{Id: uint32(len(sg.Nodes)), Name: chr, Seq: dna.RemoveBase(dna.StringToBases(v.Alt), dna.N)}
	AddNode(sg, &curr)
	return &curr
}

func isINS(v *vcf.Vcf) bool {
	var truth bool = false
	if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
		truth = true
	}
	return truth
}

func createDEL(sg *SimpleGraph, v *vcf.Vcf, chr string) *Node {
	curr := Node{Id: uint32(len(sg.Nodes)), Name: chr, Seq: dna.RemoveBase(dna.StringToBases(v.Ref[1:len(v.Ref)]), dna.N)}
	AddNode(sg, &curr)
	return &curr
}

func isDEL(v *vcf.Vcf) bool {
	var truth bool = false
	if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {
		truth = true
	}
	return truth
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
	var seeds []*SeedDev = findSeedsInMapDev(seedHash, read, seedLen, stepSize, true)
	extendSeedsDev(seeds, gg, read)
	revCompRead := fastq.Copy(read)
	fastq.ReverseComplement(revCompRead)
	var revCompSeeds []*SeedDev = findSeedsInMapDev(seedHash, revCompRead, seedLen, stepSize, false)
	extendSeedsDev(revCompSeeds, gg, revCompRead)
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
			currBest.Cigar = AddSClip(int(minQuery), len(currRead.Seq), currBest.Cigar)
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
			currBest.Cigar = AddSClip(minQuery, len(currRead.Seq), currBest.Cigar)
			currBest.Extra = PathToString(bestPath, gg)
		}
	}
	return &currBest
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


func LoopSeqToGraph(reference []*fasta.Fasta, vcfs []*vcf.Vcf) *SimpleGraph {
	gg := NewGraph()
	vcfSplit := vcf.VcfSplit(vcfs, reference)

	if len(vcfSplit) != len(reference) {
		log.Fatal("Slice of vcfs do not equal reference length")
	} else {
		for i := 0; i < len(reference); i++ {
			gg = SeqToGraph(gg, reference[i], vcfSplit[i])
		}
	}
	return gg

}

func SeqToGraph(gg *SimpleGraph, chr *fasta.Fasta, vcfFile []*vcf.Vcf) *SimpleGraph {

	var curr *Node
	var currMatch *Node
	var prev *Node
	var lastMatch *Node = nil
	var lastPos int64 = 0
	var snps []*vcf.Vcf
	for i := 0; i < len(vcfFile); i++ {
		if isSNP(vcfFile[i]) {
			snps = append(snps, vcfFile[i])
			//refAllele := &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Ref)}
			//AddNode(gg, refAllele)
			//altAllele := &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Alt)}
			//AddNode(gg, altAllele)
		}
		if isINS(vcfFile[i]) {

			if len(snps) > 0 {

				currMatch, lastPos = snpsToNodes(snps, chr, gg, lastPos)
				currMatch.Seq = append(currMatch.Seq, chr.Seq[lastPos:vcfFile[i].Pos]...)
			} else {
				currMatch = &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: chr.Seq[lastPos:vcfFile[i].Pos]}
				AddNode(gg, currMatch)
			}
			curr = &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfFile[i].Alt)}
			AddNode(gg, curr)
			lastPos++
			if lastMatch != nil {
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 0.5)
			}
			AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		if isDEL(vcfFile[i]) {
			currMatch = &Node{}
			if len(snps) > 0 {

				currMatch, lastPos = snpsToNodes(snps, chr, gg, lastPos)
				currMatch.Seq = append(currMatch.Seq, chr.Seq[lastPos:vcfFile[i].Pos]...)
			} else {
				currMatch = &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: chr.Seq[lastPos:vcfFile[i].Pos]}
				AddNode(gg, currMatch)
			}
			curr = &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfFile[i].Ref)}
			AddNode(gg, curr)
			lastPos = vcfFile[i].Pos + int64(len(curr.Seq)) - 1
			if lastMatch != nil {
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 0.5)
			}
			AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		prev = curr

	}
	lastNode := &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: chr.Seq[lastPos:len(chr.Seq)]}
	AddNode(gg, lastNode)
	AddEdge(lastMatch, lastNode, 0.5)
	AddEdge(prev, lastNode, 0.5)
	return gg
}

func snpsToNodes(snps []*vcf.Vcf, chr *fasta.Fasta, gg *SimpleGraph, lastPos int64) (*Node, int64) {
	var prevRef *Node = nil
	var prevAlt *Node = nil
	var currMatch *Node
	//refAllele := &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Ref)}
			//AddNode(gsw, refAllele)
			//altAllele := &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snp.Alt)}
			//AddNode(gsw, altAllele)
	for i := 0; i < len(snps);i++ {
		currMatch = &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: chr.Seq[lastPos:snps[i].Pos-1]}
		AddNode(gg, currMatch)
		refAllele := Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snps[i].Ref)}
		AddNode(gg, &refAllele)
		altAllele := Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snps[i].Alt)}
		AddNode(gg, &altAllele)
		AddEdge(currMatch, &refAllele, 0.5)
		AddEdge(currMatch, &altAllele, 0.5)
		if  prevRef != nil {
			AddEdge(prevRef, currMatch, 0.5)
		}
		if prevAlt != nil {
			AddEdge(prevAlt, currMatch, 0.5)
		}
		prevRef = &refAllele
		prevAlt = &altAllele
		lastPos = snps[i].Pos
	}
	currMatch = &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name}
	AddNode(gg, currMatch)
	AddEdge(prevRef, currMatch, 0.5)
	AddEdge(prevAlt, currMatch, 0.5)
	return currMatch, lastPos
}


/*
func NodesToGraph(sg *SimpleGraph, chr *fasta.Fasta, vcfFlag int, v *vcf.Vcf, idx int64, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (int64, int, *Node, *Node, *Node, *Node, *Node, *Node) {
	if chr.Name != v.Chr {
		log.Fatalf("Fasta %s does not match vcf name %s", chr.Name, v.Chr)
	}

	if vcfFlag == -1 {

		currMatch = NodeSplitByNs(sg, lastMatch, chr, idx, v.Pos)

		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {
			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)
			AddEdge(currMatch, refAllele, 0.5)
			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddEdge(currMatch, altAllele, 0.5)
			AddNode(sg, altAllele)
			vcfFlag = 1
			idx = v.Pos - 1
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			//lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			//AddNode(sg, lastMatch)
			//lastMatch = NodeSplitByNs(sg, lastMatch, chr.Seq, idx, v.Pos)

			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			prev.Seq = prev.Seq[1:]
			AddNode(sg, prev)
			AddEdge(currMatch, prev, 0.5)
			idx = v.Pos
			vcfFlag = 2
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {
			//lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			//AddNode(sg, lastMatch)
			//lastMatch = NodeSplitByNs(sg, lastMatch, chr.Seq, idx, v.Pos)
			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, prev)
			AddEdge(currMatch, prev, 0.5)

			idx = v.Pos + int64(len(prev.Seq)) - 1
			lastMatch = currMatch
			vcfFlag = 3
		}
	} else {
		if v.Pos-1-idx > 0 {
			//currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name}
			if lastMatch != nil && vcfFlag != 1 {
				AddEdge(lastMatch, currMatch, 0.5)
			}
			if vcfFlag == 1 {
				AddEdge(refAllele, currMatch, 1)
				AddEdge(altAllele, currMatch, 1)
			}
			if vcfFlag == 2 || vcfFlag == 3 {
				AddEdge(prev, currMatch, 1)
			}

			currMatch = NodeSplitByNs(sg, lastMatch, chr, idx, v.Pos)

			vcfFlag = 0
		}
		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {

			if vcfFlag == 0 {
				currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				AddEdge(currMatch, refAllele, 0.5)

				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)
				AddEdge(currMatch, altAllele, 0.5)

			} else if vcfFlag == 1 {
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, curr)
				AddEdge(refAllele, curr, 1)
				refAllele = curr

				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, curr)
				AddEdge(altAllele, curr, 1)
				altAllele = curr
			} else if vcfFlag == 2 {

				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)

				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(lastMatch, refAllele, 0.5)
				AddEdge(prev, altAllele, 1)
				curr = refAllele

			} else if vcfFlag == 3 {
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(prev, refAllele, 1)
				AddEdge(lastMatch, altAllele, 0.5)
				curr = refAllele
			} else {
				log.Fatal("Flag was not set up correctly")
			}
			lastMatch = currMatch
			vcfFlag = 1
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			curr.Seq= curr.Seq[1:]
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

				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
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
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
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
			prev = curr
			lastMatch = currMatch
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {

			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)[1:]}

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

				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
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
				currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, currMatch)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
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

			idx = v.Pos + int64(len(curr.Seq)) - 1
		}
		prev = curr
		lastMatch = currMatch

	}

	return idx, vcfFlag, curr, prev, currMatch, lastMatch, refAllele, altAllele

}*/

/*
if vcfFlag == 0 {
			currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)
			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddNode(sg, altAllele)

			AddEdge(currMatch, refAllele, 0.5)
			AddEdge(currMatch, altAllele, 0.5)
			curr = refAllele
		} else if vcfFlag == 1 {
			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, curr)
			AddEdge(refAllele, curr, 1)
			refAllele = curr

			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddNode(sg, curr)
			AddEdge(altAllele, curr, 1)
			altAllele = curr
		} else if vcfFlag == 2 {

			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)

			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddNode(sg, altAllele)

			AddEdge(lastMatch, refAllele, 0.5)
			AddEdge(prev, altAllele, 1)
			curr = refAllele

		} else if vcfFlag == 3 {
			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)
			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddNode(sg, altAllele)

			AddEdge(prev, refAllele, 1)
			AddEdge(lastMatch, altAllele, 0.5)
			curr = refAllele
		} else {
			log.Fatal("Flag was not set up correctly")
		}
*/
