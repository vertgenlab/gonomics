package genomeGraph

import (
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
)

func VariantGraph(ref <-chan fasta.Fasta, vcfMap map[string][]vcf.Vcf) *GenomeGraph {
	gg := EmptyGraph()
	var filterVcf []vcf.Vcf = make([]vcf.Vcf, 0)
	for val := range ref {
		chr := val // not sure this is necessary, but I want to make sure the function does not break
		// if the fasta being a pointer was important.
		filterVcf = vcfMap[chr.Name]
		if len(filterVcf) != 0 {
			vcf.Sort(filterVcf)
			gg = vChrGraph(gg, chr, filterVcf)
		} else {
			// Given the input is a vcf containing large structural variance
			// It is possible for a chromosome to contain no call variants (ex. chrM).
			// In such a case, we add the entire chromosome as one node.
			chrNode := &Node{Id: uint32(len(gg.Nodes)), Seq: chr.Seq, Prev: nil, Next: nil}
			AddNode(gg, chrNode)
		}
	}
	gg = SortGraph(gg)
	return gg
}

/*
func SplitGraphChr(reference []fasta.Fasta, vcfs []*vcf.Vcf) map[string]*GenomeGraph {
	gg := make(map[string]*GenomeGraph)
	vcfSplit := vcf.VcfSplit(vcfs, reference)
	if len(vcfSplit) != len(reference) {
		log.Fatal("Slice of vcfs do not equal reference length")
	} else {
		for i := 0; i < len(reference); i++ {
			gg[reference[i].Name] = vChrGraph(EmptyGraph(), reference[i], vcfSplit[i])
		}
	}
	return gg
}
*/

func vChrGraph(genome *GenomeGraph, chr fasta.Fasta, vcfsChr []vcf.Vcf) *GenomeGraph {
	vcfsChr = append(vcfsChr, vcf.Vcf{Chr: chr.Name, Pos: len(chr.Seq)})
	//log.Printf("Found %d variants on %s", len(vcfsChr), chr.Name)
	fasta.ToUpper(chr)
	var currMatch *Node = &Node{}
	var lastMatch *Node = &Node{}
	var refAllele, altAllele *Node = &Node{}, &Node{}
	var prev []*Node = nil
	var i, j, edge int
	var index int = 0
	for i = 0; i < len(vcfsChr)-1; i++ {
		if strings.Compare(chr.Name, vcfsChr[i].Chr) != 0 {
			log.Fatalf("Error: chromosome names do not match...\n")
		}
		if vcfsChr[i].Pos-index > 0 {
			currMatch = &Node{Id: uint32(len(genome.Nodes)), Seq: chr.Seq[index : vcfsChr[i].Pos-1], Prev: nil, Next: make([]Edge, 0, 2)}
			if len(currMatch.Seq) == 0 {
				currMatch = lastMatch
				//assuming we already created the ref allele, we only need to create
				//the alt alleles this iteration
				if vcf.Snp(vcfsChr[i]) {
					altAllele = &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Alt[0]), Prev: nil, Next: nil}
					altAllele = AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 0.5)
				} else if vcf.Ins(vcfsChr[i]) {
					insertion := &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Alt[0])[1:], Prev: nil, Next: nil}
					insertion = AddNode(genome, insertion)
					AddEdge(currMatch, insertion, 1)
					index = vcfsChr[i].Pos - 1
				} else if vcf.Del(vcfsChr[i]) {
					deletion := &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Ref)[1:], Prev: nil, Next: nil}
					deletion = AddNode(genome, deletion)
					AddEdge(currMatch, deletion, 1)
					if strings.Contains(vcfsChr[i].Id, "pbsv") {
						index = numbers.Min(vcfsChr[i].Pos+len(deletion.Seq)-1, vcfsChr[i+1].Pos-1)
					} else {
						index = vcfsChr[i].Pos + len(deletion.Seq)
					}
				} else if strings.Contains(vcfsChr[i].Info, "SVTYPE=SNP;INS") || strings.Contains(vcfsChr[i].Info, "SVTYPE=SNP;DEL") || strings.Contains(vcfsChr[i].Info, "SVTYPE=HAP") {
					altAllele := &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Alt[0]), Prev: nil, Next: nil}
					altAllele = AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 1)
					index = vcfsChr[i].Pos + len(refAllele.Seq) - 1
				}
				lastMatch = currMatch
			} else if len(currMatch.Seq) > 0 {
				currMatch = AddNode(genome, currMatch)
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
					refAllele = &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil}
					refAllele = AddNode(genome, refAllele)
					AddEdge(currMatch, refAllele, 0.5)

					altAllele = &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Alt[0]), Prev: nil, Next: nil}
					altAllele = AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 0.5)

					currMatch = refAllele
					index = vcfsChr[i].Pos
					for j = i + 1; j < len(vcfsChr)-1; j++ {
						if vcf.Snp(vcfsChr[j-1]) && vcf.Snp(vcfsChr[j]) && vcfsChr[j].Pos-1 == vcfsChr[j-1].Pos {
							refAllele.Seq = append(refAllele.Seq, dna.StringToBases(vcfsChr[j].Ref)...)
							altAllele.Seq = append(altAllele.Seq, dna.StringToBases(vcfsChr[j].Alt[0])...)
							index = vcfsChr[j].Pos
						} else {
							lastMatch = currMatch
							i = j - 1
							break
						}
					}
				} else if vcf.Ins(vcfsChr[i]) {
					//currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					insertion := &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Alt[0]), Prev: nil, Next: nil}
					insertion = AddNode(genome, insertion)
					AddEdge(currMatch, insertion, 1)
					index = vcfsChr[i].Pos - 1
				} else if vcf.Del(vcfsChr[i]) {
					//TODO: does not deal with ends of large deletions well, might contain extra sequence towards the end.
					//currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Alt)...)
					deletion := &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil}
					deletion = AddNode(genome, deletion)
					AddEdge(currMatch, deletion, 1)
					if strings.Contains(vcfsChr[i].Id, "pbsv") {
						index = numbers.Min(vcfsChr[i].Pos+len(deletion.Seq)-1, vcfsChr[i+1].Pos-1)
					} else {
						index = vcfsChr[i].Pos + len(deletion.Seq)
					}
				} else if isINV(vcfsChr[i]) {
					currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					inversion := mkInversionNode(genome, vcfsChr[i], chr)
					inversion = AddNode(genome, inversion)
					AddEdge(currMatch, inversion, 1)
					index = getSvEnd(vcfsChr[i])
				} else if isCNV(vcfsChr[i]) || isDup(vcfsChr[i]) {
					currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					copyVariance := &Node{Id: uint32(len(genome.Nodes)), Seq: chr.Seq[vcfsChr[i].Pos:getSvEnd(vcfsChr[i])], Prev: nil, Next: nil}
					copyVariance = AddNode(genome, copyVariance)
					AddEdge(currMatch, copyVariance, 1)
					index = getSvEnd(vcfsChr[i])
				} else if isHaplotypeBlock(vcfsChr[i]) {
					refAllele = &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Ref), Prev: nil, Next: nil}
					refAllele = AddNode(genome, refAllele)
					AddEdge(currMatch, refAllele, 1)
					altAllele = &Node{Id: uint32(len(genome.Nodes)), Seq: dna.StringToBases(vcfsChr[i].Alt[0]), Prev: nil, Next: nil}
					altAllele = AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 1)
					index = numbers.Min(vcfsChr[i].Pos+len(refAllele.Seq)-1, vcfsChr[i+1].Pos-1)
					currMatch = refAllele
				}
				lastMatch = currMatch
			}
		}
	}
	//Case: last node
	lastNode := &Node{Id: uint32(len(genome.Nodes)), Seq: chr.Seq[index:], Prev: make([]Edge, 0, len(prev)), Next: nil}
	lastNode = AddNode(genome, lastNode)
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

/*
func chrSplitByNs(chr fasta.Fasta) []fasta.Fasta {
	unGapped := bed.UngappedRegionsFromFa(chr)
	var answer []fasta.Fasta = make([]fasta.Fasta, len(unGapped))
	for i := 0; i < len(unGapped); i++ {
		answer[i] = fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d", unGapped[i].Chrom, unGapped[i].ChromStart, unGapped[i].ChromEnd), Seq: chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]}
	}
	return answer
}
*/

/*
func FaSplitByNs(fa []fasta.Fasta) []fasta.Fasta {
	var answer []fasta.Fasta
	for i := 0; i < len(fa); i++ {
		answer = append(answer, chrSplitByNs(fa[i])...)
	}
	return answer
}
*/

// TODO move these vcf helper functions to vcf
// new nodes are treated as insertion.
func isINV(v vcf.Vcf) bool {
	var truth bool = false
	data := strings.Split(v.Info, ";")
	if strings.Compare(v.Alt[0], "<INV>") == 0 || strings.Compare(data[0], "SVTYPE=INV") == 0 {
		truth = true
	}
	return truth
}

func isDup(v vcf.Vcf) bool {
	var truth bool = false
	if strings.Contains(v.Info, "SVTYPE=DUP") {
		truth = true
	}
	return truth
}

func isCNV(v vcf.Vcf) bool {
	var truth bool = false
	if strings.Contains(v.Info, "SVTYPE=CNV") {
		truth = true
	}
	return truth
}

func getSvEnd(v vcf.Vcf) int {
	if !strings.Contains(v.Info, "END=") {
		log.Fatalf("Error: Vcf might not be from PBSV...")
	} else {
		words := strings.Split(v.Info, ";")
		for i := 0; i < len(words); i++ {
			if strings.Contains(words[i], "END=") {
				text := strings.Split(words[i], "END=")
				return common.StringToInt(text[1])
			}
		}
	}
	return 0
}

func mkInversionNode(genome *GenomeGraph, v vcf.Vcf, chr fasta.Fasta) *Node {
	invSeq := make([]dna.Base, 0)
	invSeq = append(invSeq, chr.Seq[v.Pos:getSvEnd(v)]...)
	dna.ReverseComplement(invSeq)
	inversion := &Node{Id: uint32(len(genome.Nodes)), Seq: invSeq, Prev: nil, Next: nil}
	return inversion
}

/*
func createSNP(sg *GenomeGraph, snp *vcf.Vcf, chr string) (*Node, *Node) {
	refAllele := &Node{Id: uint32(len(sg.Nodes)), Seq: dna.StringToBases(snp.Ref), Next: nil, Prev: nil}
	altAllele := &Node{Id: uint32(len(sg.Nodes) + 1), Seq: dna.StringToBases(snp.Alt[0]), Next: nil, Prev: nil}
	AddNode(sg, refAllele)
	AddNode(sg, altAllele)
	return refAllele, altAllele
}
*/

/*
func NodeSplitByNs(sg *GenomeGraph, currMatch *Node, chr *fasta.Fasta, index int64, end int64) *Node {
	var inRegion bool = false
	var start int64 = 0
	for ; index < end; index++ {
		if dna.DefineBase(chr.Seq[start]) && inRegion == false {
			inRegion = false
			currMatch.Seq = chr.Seq[start:index]
			start = index
		} else if dna.DefineBase(chr.Seq[start]) && inRegion == true {
			newMatch := &Node{Id: uint32(len(sg.Nodes))}
			AddNode(sg, newMatch)
			AddEdge(currMatch, newMatch, 1)
			inRegion = true
			currMatch = newMatch
		}
	}
	return currMatch
}
*/

/*
type noNsBed struct {
	Start int32
	End   int32
}
*/

/*
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
*/

/*
func createINS(sg *GenomeGraph, v *vcf.Vcf, chr string) *Node {
	curr := Node{Id: uint32(len(sg.Nodes)), Seq: dna.StringToBases(v.Alt[0])}
	AddNode(sg, &curr)
	return &curr
}
*/

func isHaplotypeBlock(v vcf.Vcf) bool {
	if strings.Contains(v.Info, "SVTYPE=SNP;INS") || strings.Contains(v.Info, "SVTYPE=SNP;DEL") || strings.Contains(v.Info, "SVTYPE=HAP") {
		return true
	} else {
		return false
	}
}

/*
func createDEL(sg *GenomeGraph, v *vcf.Vcf, chr string) *Node {
	curr := Node{Id: uint32(len(sg.Nodes)), Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}
	AddNode(sg, &curr)
	return &curr
}
*/
