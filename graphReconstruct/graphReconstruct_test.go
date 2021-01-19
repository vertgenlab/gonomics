package graphReconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"testing"
)

var GCcontent = 0.42
var input = []struct {
	newickFilename string // first input
	length         int    // second input
}{
	{"testdata/newickLongBranches.txt", 1005},
}

func TestGraphRecon(t *testing.T) {
	var humanNode1 = &simpleGraph.Node{Id: 0, Name: "humanNode1", Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	var humanNode2 = &simpleGraph.Node{Id: 1, Name: "humanNode2", Seq: dna.StringToBases("AAA"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("AAA")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	var humanNode3 = &simpleGraph.Node{Id: 2, Name: "humanNode3", Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	var humanNode4 = &simpleGraph.Node{Id: 3, Name: "humanNode4", Seq: dna.StringToBases("CCC"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("CCC")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	var humanNode5 = &simpleGraph.Node{Id: 4, Name: "humanNode5", Seq: dna.StringToBases("GGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("GGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}

	var humanEdge1 = &simpleGraph.Edge{humanNode2, 0.25}
	var humanEdge2 = &simpleGraph.Edge{humanNode3, 0.75}
	var humanEdge3 = &simpleGraph.Edge{humanNode3, 1.00}
	var humanEdge4 = &simpleGraph.Edge{humanNode4, 0.25}
	var humanEdge5 = &simpleGraph.Edge{humanNode5, 0.75}
	var humanEdge6 = &simpleGraph.Edge{humanNode5, 1.00}

	humanNode1.Next = []*simpleGraph.Edge{humanEdge1, humanEdge2}
	humanNode2.Next = []*simpleGraph.Edge{humanEdge3}
	humanNode3.Next = []*simpleGraph.Edge{humanEdge4, humanEdge5}
	humanNode4.Next = []*simpleGraph.Edge{humanEdge6}

	var humanGraph = &simpleGraph.SimpleGraph{Nodes: []*simpleGraph.Node{humanNode1, humanNode2, humanNode3, humanNode4, humanNode5}}

	simpleGraph.PrintGraph(humanGraph)

	var chimpNode1 = &simpleGraph.Node{Id: 0, Name: "chimpNode1", Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	var chimpNode2 = &simpleGraph.Node{Id: 1, Name: "chimpNode2", Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}

	var chimpEdge1 = &simpleGraph.Edge{chimpNode2, 1.00}

	chimpNode1.Next = []*simpleGraph.Edge{chimpEdge1}

	var chimpGraph = &simpleGraph.SimpleGraph{Nodes: []*simpleGraph.Node{chimpNode1, chimpNode2}}

	simpleGraph.PrintGraph(chimpGraph)

	var gorillaNode1 = &simpleGraph.Node{Id: 0, Name: "gorillaNode1", Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	var gorillaNode2 = &simpleGraph.Node{Id: 1, Name: "gorillaNode2", Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}

	var gorillaEdge1 = &simpleGraph.Edge{Dest: gorillaNode2, Prob: 1.00}

	gorillaNode1.Next = []*simpleGraph.Edge{gorillaEdge1}

	var gorillaGraph = &simpleGraph.SimpleGraph{Nodes: []*simpleGraph.Node{gorillaNode1, gorillaNode2}}

	simpleGraph.PrintGraph(gorillaGraph)

	var nodeAlign0 = graphColumn{AlignId: 0, AlignNodes: [][]*simpleGraph.Node{[]*simpleGraph.Node{humanGraph.Nodes[0]}, []*simpleGraph.Node{chimpGraph.Nodes[0]}, []*simpleGraph.Node{gorillaGraph.Nodes[0]}}}
	var nodeAlign1 = graphColumn{AlignId: 1, AlignNodes: [][]*simpleGraph.Node{[]*simpleGraph.Node{humanGraph.Nodes[1]}}}
	var nodeAlign2 = graphColumn{AlignId: 2, AlignNodes: [][]*simpleGraph.Node{[]*simpleGraph.Node{humanGraph.Nodes[2]}, []*simpleGraph.Node{chimpGraph.Nodes[1]}, []*simpleGraph.Node{gorillaGraph.Nodes[1]}}}

	log.Print(nodeAlign0.AlignId)
	log.Print(nodeAlign1.AlignId)
	log.Print(nodeAlign2.AlignId)

	path, prob := PathFinder(humanGraph)
	log.Print("path")
	log.Print(path)
	log.Print("prob")
	log.Print(prob)

}

//func Test_reconstruct(t *testing.T) {
//	for _, test := range input {
//		tre, er := expandedTree.ReadNewick(test.newickFilename)
//		if er != nil {
//			log.Fatal("Couldn't read file")
//		}
//		fasta.Write("RandGeneOutput.fasta", simulate.RandGene("test", test.length, GCcontent)) //galGal6 GC
//		simulate.Simulate("RandGeneOutput.fasta", tre, "testdata/genePred.gp")
//		WriteTreeToFasta(tre, "simOut.fasta")
//		WriteLeavesToFasta(tre, "leavesOnly.Fasta")
//
//		tr := expandedTree.ReadTree(test.newickFilename, "leavesOnly.fasta")
//		leaves := expandedTree.GetLeaves(tr)
//		for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
//			LoopNodes(tr, i)
//		}
//		WriteTreeToFasta(tr, "reconOut.fasta")
//		accuracyData := ReconAccuracy("simOut.fasta", "reconOut.fasta")
//		for name, accuracy := range accuracyData {
//			log.Printf("%s %f \n", name, accuracy)
//		}
//	}
//}
