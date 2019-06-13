package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/graph"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	axtAlign := axt.Read(flag.Arg(0))
	vcfAxt := axt.CallSnpsToVcf(axtAlign)
	ref := fasta.Read(flag.Arg(1))
	gGraph := graph.NewGraph()
	g := graph.RefernceToGraph(vcfAxt, ref, gGraph)

	fmt.Println(g.String())
	vcf.PrintVcf(vcfAxt)
	//graph.PrintGraph(g)
	graph.Write("dev.gg", g)
	gg := graph.Read("dev.gg")
	graph.PrintGraph(gg)
}

func usage() {
	fmt.Print(
		"GenomeGraph dev - b\n")
	flag.PrintDefaults()
}

/*

//Modified from align package
func AlignToGraph(alpha *fasta.Fasta, beta *fasta.Fasta) *GenomeGraph {
	g := NewGraph()
	nodeCount := 0
	var curr *Node
	var prev *Node
	_, route := align.AffineGap(alpha.Seq, beta.Seq, align.DefaultScoreMatrix, -400, -30)
	var pairSeq []*fasta.Fasta
	pairSeq = append(pairSeq, alpha)
	pairSeq = append(pairSeq, beta)
	mAlign := align.AllSeqAffine(pairSeq, align.DefaultScoreMatrix, -400, -30)
	var seqIdx int64 = 0
	for i, _ := range route {
		oper := route[i].Op
		//saved := &Node{qSeq: nil}
		switch oper {
		//case align.colM:
		case 0:
			curr = &Node{qSeq: qDna.FromDnaToQFrag(mAlign[0].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[0].Name)}
			g.AddNode(curr)
			nodeCount++
			if nodeCount > 0 {
				g.AddEdge(prev, curr, 1)
			}
		//case align.colI:
		case 1:
			curr = &Node{qSeq: qDna.FromDnaToQFrag(mAlign[1].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[1].Name)}
			g.AddNode(curr)
			nodeCount++
			if nodeCount > 1 {
				g.AddEdge(prev, curr, 1)
			}
		//case align.colD:
		case 2:
			curr = &Node{qSeq: qDna.FromDnaToQFrag(mAlign[0].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[0].Name)}
			g.AddNode(curr)
			nodeCount++
			if nodeCount > 0 {
				g.AddEdge(prev, curr, 1)
			}
		}
		seqIdx = seqIdx + route[i].RunLength
		prev = curr
	}
	fmt.Println(route)
	fmt.Println(dna.BasesToString(mAlign[0].Seq))
	fmt.Println(dna.BasesToString(mAlign[1].Seq))
	return g
}

func newGenomeGraph(inFile string) {
	records, err := fasta.ReadNew(inFile)
	for j, _ := range records {
		fmt.Println(dna.BasesToString(records[j].Seq))
	}
	a := fasta.DivideFastaAll(records, 3)
	g := FillGraph(a)
	if err != nil {
		log.Fatal(err)
	}
	s := g.String()
	fmt.Println(s)
	//mAlign := align.AllSeqAffine(records, align.DefaultScoreMatrix, -400, -30)
	alpha := records[0]
	beta := records[1]
	//z, cigar := align.AffineGap(alpha, beta, align.DefaultScoreMatrix, -400, -30)
	//fmt.Println(dna.BasesToString(mAlign[0].Seq))
	//fmt.Println(dna.BasesToString(mAlign[1].Seq))
	//fmt.Println(z)
	//fmt.Println(cigar)
	gg := AlignToGraph(alpha, beta)
	fmt.Println(gg.String())

}
*/

/*
func FillGraph(input [][]*fasta.Fasta) *GenomeGraph {
	g := NewGraph()
	var prev *Node
	var curr *Node
	for i := 0; i < len(input); i++ {
		for j := 0; j < len(input[i]); j++ {
			curr = &Node{qSeq: qDna.FromFasta(input[i][j])}
			g.AddNode(curr)
			if j > 0 {
				g.AddEdge(prev, curr, 1)
			}
			prev = curr
		}
	}
	return g
}
*/
