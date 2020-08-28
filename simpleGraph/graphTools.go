package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

func VariantGraph(ref <-chan *fasta.Fasta, vcfMap map[string][]*vcf.Vcf) *SimpleGraph {
	gg := NewGraph()
	var filterVcf []*vcf.Vcf = make([]*vcf.Vcf, 0)
	for chr := range ref {
		filterVcf = vcfMap[chr.Name]
		if len(filterVcf) != 0 {
			vcf.Sort(filterVcf)
			gg = vChrGraph(gg, chr, filterVcf)
		} else {
			// Given the input is a vcf containing large structural variance
			// It is possible for a chromosome to contain no call variants (ex. chrM).
			// In such a case, we add the entire chromosome as one node.
			chrNode := &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, Seq: chr.Seq, Prev: nil, Next: nil, Info: Annotation{Allele: 0, Start: 1, Variant: 0}}
			AddNode(gg, chrNode)
		}
	}
	gg = SortGraph(gg)
	return gg
}

func SplitGraphChr(reference []*fasta.Fasta, vcfs []*vcf.Vcf) map[string]*SimpleGraph {
	gg := make(map[string]*SimpleGraph)
	vcfSplit := vcf.VcfSplit(vcfs, reference)
	if len(vcfSplit) != len(reference) {
		log.Fatal("Slice of vcfs do not equal reference length")
	} else {
		for i := 0; i < len(reference); i++ {
			gg[reference[i].Name] = vChrGraph(NewGraph(), reference[i], vcfSplit[i])
		}
	}
	return gg
}

func vChrGraph(genome *SimpleGraph, chr *fasta.Fasta, vcfsChr []*vcf.Vcf) *SimpleGraph {
	vcfsChr = append(vcfsChr, &vcf.Vcf{Chr: chr.Name, Pos: int64(len(chr.Seq))})
	//log.Printf("Found %d variants on %s", len(vcfsChr), chr.Name)
	fasta.ToUpper(chr)
	var currMatch *Node = &Node{}
	var lastMatch *Node = &Node{}
	var refAllele, altAllele *Node = &Node{}, &Node{}
	var prev []*Node = nil
	var i, j, edge int
	var index int64 = 0
	for i = 0; i < len(vcfsChr)-1; i++ {
		if strings.Compare(chr.Name, vcfsChr[i].Chr) != 0 {
			log.Fatalf("Error: chromosome names do not match...\n")
		}
		if vcfsChr[i].Pos-index > 0 {
			currMatch = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: chr.Seq[index : vcfsChr[i].Pos-1], Prev: nil, Next: make([]*Edge, 0, 2), Info: Annotation{Allele: 0, Start: uint32(index + 1), Variant: 0}}
			if len(currMatch.Seq) == 0 {
				currMatch = lastMatch
				//assuming we already created the ref allele, we only need to create
				//the alt alleles this iteration
				if vcf.Snp(vcfsChr[i]) {
					altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Alt), Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 1}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 0.5)
				} else if vcf.Ins(vcfsChr[i]) {
					insertion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Alt)[1:], Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 2}}
					AddNode(genome, insertion)
					AddEdge(currMatch, insertion, 1)

					index = vcfsChr[i].Pos - 1
				} else if vcf.Del(vcfsChr[i]) {
					deletion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Ref)[1:], Prev: nil, Next: nil, Info: Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 3}}
					AddNode(genome, deletion)
					AddEdge(currMatch, deletion, 1)
					if strings.Contains(vcfsChr[i].Id, "pbsv") {
						index = common.MinInt64(vcfsChr[i].Pos+int64(len(deletion.Seq))-1, vcfsChr[i+1].Pos-1)
					} else {
						index = vcfsChr[i].Pos + int64(len(deletion.Seq))
					}
				} else if strings.Compare(vcfsChr[i].Format, "SVTYPE=SNP;INS") == 0 || strings.Compare(vcfsChr[i].Format, "SVTYPE=SNP;DEL") == 0 || strings.Compare(vcfsChr[i].Format, "SVTYPE=HAP") == 0 {
					altAllele := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Alt), Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 4}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 1)
					index = vcfsChr[i].Pos + int64(len(refAllele.Seq)) - 1
				}
				lastMatch = currMatch
			} else if len(currMatch.Seq) > 0 {

				AddNode(genome, currMatch)
				if lastMatch != nil {
					if len(lastMatch.Next) > 0 {
						for edge = 0; edge < len(lastMatch.Next); edge++ {
							AddEdge(lastMatch.Next[edge].Dest, currMatch, 1)
						}
					}
					if i > 0 {
						if vcf.Snp(vcfsChr[i-1]) || isHaplotypeBlock(vcfsChr[i-1]) {
							AddEdge(altAllele, currMatch, 1)
						}
					}
					AddEdge(lastMatch, currMatch, 1)
					SetEvenWeights(lastMatch)
				}
				if vcf.Snp(vcfsChr[i]) {

					refAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil, Info: Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 1}}
					AddNode(genome, refAllele)
					AddEdge(currMatch, refAllele, 0.5)

					altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Alt), Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 1}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 0.5)

					currMatch = refAllele
					index = vcfsChr[i].Pos
					for j = i + 1; j < len(vcfsChr)-1; j++ {
						if vcf.Snp(vcfsChr[j-1]) && vcf.Snp(vcfsChr[j]) && vcfsChr[j].Pos-1 == vcfsChr[j-1].Pos {
							refAllele.Seq = append(refAllele.Seq, dna.StringToBases(vcfsChr[j].Ref)...)
							altAllele.Seq = append(altAllele.Seq, dna.StringToBases(vcfsChr[j].Alt)...)
							index = vcfsChr[j].Pos
						} else {
							lastMatch = currMatch
							i = j - 1
							break
						}
					}
				} else if vcf.Ins(vcfsChr[i]) {
					//currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					insertion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Alt), Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 2}}
					AddNode(genome, insertion)
					AddEdge(currMatch, insertion, 1)
					index = vcfsChr[i].Pos - 1
				} else if vcf.Del(vcfsChr[i]) {
					//TODO: does not deal with ends of large deletions well, might contain extra sequence towards the end.
					//currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Alt)...)
					deletion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil, Info: Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 3}}
					AddNode(genome, deletion)
					AddEdge(currMatch, deletion, 1)
					if strings.Contains(vcfsChr[i].Id, "pbsv") {
						index = common.MinInt64(vcfsChr[i].Pos+int64(len(deletion.Seq))-1, vcfsChr[i+1].Pos-1)
					} else {
						index = vcfsChr[i].Pos + int64(len(deletion.Seq))
					}
				} else if isINV(vcfsChr[i]) {
					currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					inversion := mkInversionNode(genome, vcfsChr[i], chr)
					AddNode(genome, inversion)
					AddEdge(currMatch, inversion, 1)
					index = getSvEnd(vcfsChr[i])
				} else if isCNV(vcfsChr[i]) || isDup(vcfsChr[i]) {
					currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					copyVariance := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: chr.Seq[vcfsChr[i].Pos:getSvEnd(vcfsChr[i])], Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 6}}
					AddNode(genome, copyVariance)
					AddEdge(currMatch, copyVariance, 1)
					index = getSvEnd(vcfsChr[i])
				} else if isHaplotypeBlock(vcfsChr[i]) {
					refAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil, Info: Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 4}}
					AddNode(genome, refAllele)
					AddEdge(currMatch, refAllele, 1)
					altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfsChr[i].Alt), Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 4}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 1)
					index = common.MinInt64(vcfsChr[i].Pos+int64(len(refAllele.Seq))-1, vcfsChr[i+1].Pos-1)
					currMatch = refAllele
				}
				lastMatch = currMatch
			} else {
				//num++
				//log.Printf("index=%d, vcfPos=%d\ndiff=%d\n%s", index, vcfsChr[i].Pos, vcfsChr[i].Pos-index, vcfsChr[i].Info)
			}
		}
	}
	//Case: last node
	lastNode := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: chr.Seq[index:], Prev: make([]*Edge, 0, len(prev)), Next: nil, Info: Annotation{Allele: 0, Start: uint32(index + 1), Variant: 0}}
	AddNode(genome, lastNode)
	for edge = 0; edge < len(lastMatch.Next); edge++ {
		AddEdge(lastMatch.Next[edge].Dest, lastNode, 1)
	}
	if vcf.Snp(vcfsChr[len(vcfsChr)-2]) || isHaplotypeBlock(vcfsChr[len(vcfsChr)-2]) {
		AddEdge(altAllele, lastNode, 1)
	}
	AddEdge(lastMatch, lastNode, 1)
	SetEvenWeights(lastMatch)
	return genome
}

func GraphToFa(gg *SimpleGraph) []*fasta.Fasta {
	var answer []*fasta.Fasta
	for i := 0; i < len(gg.Nodes); i++ {
		if len(gg.Nodes[i].Seq) == 0 {
			log.Printf("Empty node: %s_%d_%d_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start)
		} else {
			fragment := &fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start), Seq: gg.Nodes[i].Seq}
			answer = append(answer, fragment)
		}
	}
	return answer
}

func chrSplitByNs(chr *fasta.Fasta) []*fasta.Fasta {
	unGapped := bed.UngappedRegionsFromFa(chr)
	var answer []*fasta.Fasta = make([]*fasta.Fasta, len(unGapped))
	for i := 0; i < len(unGapped); i++ {
		answer[i] = &fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d", unGapped[i].Chrom, unGapped[i].ChromStart, unGapped[i].ChromEnd), Seq: chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]}
	}
	return answer
}

func FaSplitByNs(fa []*fasta.Fasta) []*fasta.Fasta {
	var answer []*fasta.Fasta
	for i := 0; i < len(fa); i++ {
		answer = append(answer, chrSplitByNs(fa[i])...)
	}
	return answer
}

/*
func ToFaSubSet(gg *SimpleGraph, variant uint8) []*fasta.Fasta {
	var answer []*fasta.Fasta
	for i := 0; i < len(gg.Nodes); i++ {
		if len(gg.Nodes[i].Seq) == 0 {
			log.Printf("Empty node: %s_%d_%d_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start)
		} else {
			if gg.Nodes[i].Info.Variant == variant {
				fragment := &fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start), Seq: gg.Nodes[i].Seq}
				answer = append(answer, fragment)
			}
		}
	}
	return answer
}*/
//TODO move these vcf helper functions to vcf
//new nodes are treated as insertion
func isINV(v *vcf.Vcf) bool {
	var truth bool = false
	data := strings.Split(v.Info, ";")
	if strings.Compare(v.Alt, "<INV>") == 0 || strings.Compare(data[0], "SVTYPE=INV") == 0 {
		truth = true
	}
	return truth
}

func isDup(v *vcf.Vcf) bool {
	var truth bool = false
	if strings.Contains(v.Info, "SVTYPE=DUP") {
		truth = true
	}
	return truth
}

func isCNV(v *vcf.Vcf) bool {
	var truth bool = false
	if strings.Contains(v.Info, "SVTYPE=CNV") {
		truth = true
	}
	return truth
}

func getSvEnd(v *vcf.Vcf) int64 {
	if !strings.Contains(v.Info, "END=") {
		log.Fatalf("Error: Vcf might not be from PBSV...")
	} else {
		words := strings.Split(v.Info, ";")
		for i := 0; i < len(words); i++ {
			if strings.Contains(words[i], "END=") {
				text := strings.Split(words[i], "END=")
				return common.StringToInt64(text[1])
			}
		}
	}
	return 0
}

func mkInversionNode(genome *SimpleGraph, v *vcf.Vcf, chr *fasta.Fasta) *Node {
	invSeq := make([]dna.Base, 0)
	invSeq = append(invSeq, chr.Seq[v.Pos:getSvEnd(v)]...)
	dna.ReverseComplement(invSeq)
	inversion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, Seq: invSeq, Prev: nil, Next: nil, Info: Annotation{Allele: 1, Start: uint32(v.Pos), Variant: 5}}
	return inversion
}

func createSNP(sg *SimpleGraph, snp *vcf.Vcf, chr string) (*Node, *Node) {
	refAllele := &Node{Id: uint32(len(sg.Nodes)), Name: chr, Seq: dna.StringToBases(snp.Ref), Next: nil, Prev: nil}
	altAllele := &Node{Id: uint32(len(sg.Nodes) + 1), Name: chr, Seq: dna.StringToBases(snp.Alt), Next: nil, Prev: nil}
	AddNode(sg, refAllele)
	AddNode(sg, altAllele)
	return refAllele, altAllele
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

func createINS(sg *SimpleGraph, v *vcf.Vcf, chr string) *Node {
	curr := Node{Id: uint32(len(sg.Nodes)), Name: chr, Seq: dna.StringToBases(v.Alt)}
	AddNode(sg, &curr)
	return &curr
}

func isHaplotypeBlock(v *vcf.Vcf) bool {
	if strings.Compare(v.Format, "SVTYPE=SNP;INS") == 0 || strings.Compare(v.Format, "SVTYPE=SNP;DEL") == 0 || strings.Compare(v.Format, "SVTYPE=HAP") == 0 {
		return true
	} else {
		return false
	}
}

func createDEL(sg *SimpleGraph, v *vcf.Vcf, chr string) *Node {
	curr := Node{Id: uint32(len(sg.Nodes)), Name: chr, Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}
	AddNode(sg, &curr)
	return &curr
}

/*
func FaToGenomeGraph(ref []*fasta.Fasta, vcfs []*vcf.Vcf) (*SimpleGraph, map[string]*Node) {
	gg := NewGraph()
	chrom := make(map[string]*Node)

	vcfSplit := vcf.VcfSplit(vcfs, ref)
	var i, nodeIdx int
	if len(vcfSplit) != len(ref) {
		log.Fatal("Slice of vcfs do not equal reference length")
	} else {
		for i = 0; i < len(ref); i++ {
			gg = vChrGraph(gg, ref[i], vcfSplit[i])
			_, ok := chrom[ref[i].Name]
			if !ok {
				chrom[ref[i].Name] = gg.Nodes[nodeIdx]
				nodeIdx = len(gg.Nodes)
			} else {
				log.Fatalf("Error: The indexing to keep track of the head nodes are off...")
			}
		}
	}
	return gg, chrom
}*/
