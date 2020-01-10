package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

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

func NodesToGraph(sg *SimpleGraph, chr *fasta.Fasta, vcfFlag int, v *vcf.Vcf, idx int64, curr *Node, prev *Node, currMatch *Node, lastMatch *Node, refAllele *Node, altAllele *Node) (int64, int, *Node, *Node, *Node, *Node, *Node, *Node) {
	if vcfFlag == -1 && v.Pos == 1 {
		if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {
			refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, refAllele)
			altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}

			AddNode(sg, altAllele)
			vcfFlag = 1
			idx++
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
			AddNode(sg, lastMatch)
			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}
			AddNode(sg, prev)
			AddEdge(lastMatch, prev, 0.5)
			idx = v.Pos
			vcfFlag = 2
		}
		if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {
			lastMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
			AddNode(sg, lastMatch)
			prev = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt[1:len(v.Ref)])}
			AddNode(sg, prev)
			AddEdge(lastMatch, prev, 0.5)

			idx = v.Pos + int64(len(prev.Seq)) - 1
			vcfFlag = 3
		}
	} else {
		if v.Pos-1-idx > 0 {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:v.Pos]}
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
			currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
			if vcfFlag == 0 {
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref)}
				AddNode(sg, refAllele)
				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt)}
				AddNode(sg, altAllele)

				AddEdge(currMatch, refAllele, 0.5)
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
			vcfFlag = 1
			idx = v.Pos
		}
		if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}

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
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}
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
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}
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

			curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}

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
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Ref[1:len(v.Ref)])}
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
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(v.Alt[1:len(v.Alt)])}
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

/*
func SeqToGraph(vcfFile []*vcf.Vcf, sequence *fasta.Fasta, gsw *SimpleGraph) *SimpleGraph {
	g := gsw


	var curr []*Node
	var currMatch []*Node
	var prev []*Node
	var lastMatch []*Node = []*Node{}
	var lastPos int64
	//var snps []*vcf.Vcf
	lastPos = 0
	for i := 0; i < len(vcfFile); i++ {
		if strings.Compare(vcfFile[i].Format, "SVTYPE=SNP") == 0 {
			//snps = append(snps, vcfFile[i])
			curr = &Node{Id: uint32(len(g.Nodes)), sequence.Name, Seq: []dna.Base{1}}
			AddNode(g, curr)
		}
		//case insertion
		if strings.Compare(vcfFile[i].Format, "SVTYPE=INS") == 0 {
			//logic for calling SNPs

			//added location of fragment and converted to 1 base by add 1
			currMatch = &Node{QSeq: SnpQFrag(sequence.Seq[lastPos:vcfFile[i].Pos], snpSeq[lastPos:vcfFile[i].Pos], sequence.Name, lastPos+1, vcfFile[i].Pos), NodeID: g.NumNode}
			snps = nil
			g.AddNode(currMatch)
			bases := dna.StringToBases(vcfFile[i].Alt)
			curr = &Node{QSeq: qDna.QFragCoord(bases[1:len(bases)], sequence.Name, lastPos+1, lastPos+1), NodeID: g.NumNode}
			g.AddNode(curr)
			lastPos++
			if len(lastMatch) != 0 {
				g.AddEdge(lastMatch, currMatch, 1)
				g.AddEdge(prev, currMatch, 1)
			}
			g.AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		//deletion in vcf record
		if strings.Compare(vcfFile[i].Format, "SVTYPE=DEL") == 0 {


			currMatch = &Node{QSeq: SnpQFrag(sequence.Seq[lastPos:vcfFile[i].Pos], snpSeq[lastPos:vcfFile[i].Pos], sequence.Name, lastPos+1, vcfFile[i].Pos), NodeID: g.NumNode}
			snps = nil
			g.AddNode(currMatch)

			bases := dna.StringToBases(vcfFile[i].Ref)
			curr = &Node{QSeq: qDna.QFragCoord(bases[1:len(bases)], sequence.Name, vcfFile[i].Pos+1, vcfFile[i].Pos+int64(len(bases))-1), NodeID: g.NumNode}
			g.AddNode(curr)
			lastPos = vcfFile[i].Pos + int64(len(bases)) - 1
			if len(lastMatch) != 0 {
				g.AddEdge(lastMatch, currMatch, 1)
				g.AddEdge(prev, currMatch, 1)
			}
			g.AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		prev = curr
	}
	//last match case:
	snpSeq := SnpToSeq(sequence.Seq, snps)
	lastNode := &Node{QSeq: SnpQFrag(sequence.Seq[lastPos:len(sequence.Seq)], snpSeq[lastPos:len(sequence.Seq)], sequence.Name, lastPos+1, int64(len(sequence.Seq))), NodeID: g.NumNode}
	g.AddNode(lastNode)
	g.AddEdge(prev, lastNode, 1)
	g.AddEdge(lastMatch, lastNode, 1)
	return g
}

func snpToNodes(sg *SimpleGraph, chr *fasta.Fasta, snps []*vcf.Vcf) {
	var idx int64 = 0

	var match Node = Node{Id: 0, Name: "", Seq: nil}
	var refAllele Node = Node{Id: 0, Name: "", Seq: nil}
	var altAllele Node = Node{Id: 0, Name: "", Seq: nil}

	for i := 0; i < len(snps);i++ {
		if refAllele.Seq != nil && altAllele.Seq != nil {
			//make new node for ref snp
			AddNode(sg, &refAllele)

			AddNode(sg, &altAllele)

			//match = Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx+1:snps[i].Pos-2]}
			match.Seq = chr.Seq[idx+1:snps[i].Pos-2]
			AddNode(sg, &match)

			AddEdge(&refAllele, &match, 0.5)
			AddEdge(&altAllele, &match, 0.5)

			refAllele = Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snps[i].Ref)}

			altAllele = Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(snps[i].Alt)}


			idx = snps[i].Pos-1
		}

	}
	if idx < int64(len(chr.Seq)-1) {
		match.Seq = chr.Seq[idx:]
		AddEdge(&refAllele, &match, 0.5)
		AddEdge(&altAllele, &match, 0.5)
	}
	//connect last snp to last match

}

func vcfToNodes(sg *SimpleGraph, chr *fasta.Fasta, vcfs []*vcf.Vcf) *SimpleGraph {
	//g := sg
	var idx int64 = 0

	var curr *Node
	var prev *Node

	var currMatch *Node
	var lastMatch *Node = nil

	var refAllele *Node
	var altAllele *Node

	var delSeq []dna.Base

	var vcfFlag int = -1
	for i := 0; i < len(vcfs); i++ {
		if vcfFlag == -1 && vcfs[i].Pos-idx-1 > 0 {
			//no nodes were created yet, creating first matches
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:vcfs[i].Pos-2]}
			AddNode(sg, currMatch)

			vcfFlag = 0
		} else if vcfFlag == 1 {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:vcfs[i].Pos-2]}
			AddNode(sg, currMatch)

			AddEdge(refAllele, currMatch, 0.5)
			AddEdge(altAllele, currMatch, 0.5)
			vcfFlag = 0
		} else if vcfFlag == 2 || vcfFlag == 3 {
			currMatch = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:vcfs[i].Pos-2]}
			AddNode(sg, currMatch)
			if lastMatch == nil {
				AddEdge(currMatch, curr, 0.5)
			} else {
				AddEdge(lastMatch, currMatch, 0.5)
				AddEdge(prev, currMatch, 0.5)
			}
			vcfFlag = 0
		} else {
			if strings.Compare(vcfs[i].Format, "SVTYPE=SNP") == 0 {
				refAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfs[i].Ref)}
				altAllele = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfs[i].Alt)}
				if vcfFlag == 0 {


					//node operations
					AddNode(sg, refAllele)
					AddNode(sg, altAllele)
					AddEdge(currMatch, refAllele, 0.5)
					AddEdge(currMatch, altAllele, 0.5)
				} else if vcfFlag == 1 {
					//connect the consecutive SNPs
					curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfs[i].Ref)}
					AddNode(sg, curr)
					AddEdge(refAllele, curr, 1)
					refAllele = curr

					curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfs[i].Alt)}
					AddNode(sg, curr)
					AddEdge(altAllele, curr, 1)
					altAllele = curr

				} else if vcfFlag == 2 {

					AddNode(sg, refAllele)
					AddNode(sg, altAllele)

					AddEdge(currMatch, refAllele, 0.5)
					AddEdge(prev, altAllele, 1)

				} else if vcfFlag == 3 {
					AddNode(sg, refAllele)
					AddNode(sg, altAllele)

					AddEdge(prev, refAllele, 1)
					AddEdge(currMatch, altAllele, 0.5)
				} else {
					log.Fatal("Flag was not set up correctly")
				}
				vcfFlag = 1
				idx = vcfs[i].Pos+1
			}
			if strings.Compare(vcfs[i].Format, "SVTYPE=INS") == 0 {
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: dna.StringToBases(vcfs[i].Alt)}
				AddNode(sg, curr)
				if vcfFlag == 0 {
					AddEdge(currMatch, curr, 0.5)
					vcfFlag = 2
				} else if vcfFlag == 1 {
					AddEdge(altAllele, curr, 0.5)
					altAllele = curr
					//flag as SNP because there are two nodes that need to be connected to the next match
					vcfFlag = 1
				} else if vcfFlag == 2 {
					AddEdge(prev, curr, 1)
					vcfFlag = 2
				} else if vcfFlag == 3 {
					//refAllele = prev
					//altAllele = curr
					AddEdge(currMatch, curr, 0.5)
					vcfFlag = 2
				} else {
					log.Fatal("Flag was not set up correctly")
				}
				prev = curr

				idx++
			}
			if strings.Compare(vcfs[i].Format, "SVTYPE=DEL") == 0 {
				delSeq = dna.StringToBases(vcfs[i].Ref)
				curr = &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: delSeq}
				AddNode(sg, curr)
				if vcfFlag == 0 {
					AddEdge(currMatch, curr, 0.5)
				} else if vcfFlag == 1 {
					AddEdge(refAllele, curr, 0.5)
					refAllele = curr
				} else if vcfFlag == 2 {
					AddEdge(currMatch, curr, 0.5)
				} else if vcfFlag == 3 {
					AddEdge(prev, curr, 1)
				} else {
					log.Fatal("Flag was not set up correctly")
				}
				prev = curr
				idx = vcfs[i].Pos + int64(len(delSeq)) - 1
				vcfFlag = 3
			}
			lastMatch = currMatch
		}


	}
	//last node case
	lastNode := &Node{Id: uint32(len(sg.Nodes)), Name: chr.Name, Seq: chr.Seq[idx:]}
	AddNode(sg, lastNode)
	if vcfFlag == 1 {
		AddEdge(refAllele, lastNode, 0.5)
		AddEdge(altAllele, lastNode, 0.5)
	}
	if vcfFlag == 2 || vcfFlag == 3 {
		AddEdge(lastMatch, lastNode, 0.5)
		AddEdge(prev, lastNode, 0.5)
	}
	return sg

}



*/
