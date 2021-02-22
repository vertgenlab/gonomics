package graphReconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"testing"
)

var (
	allAlign []graphColumn

	humanNode1 = &simpleGraph.Node{Id: 0, Name: "humanNode1", Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	humanNode2 = &simpleGraph.Node{Id: 1, Name: "humanNode2", Seq: dna.StringToBases("AAA"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("AAA")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	humanNode3 = &simpleGraph.Node{Id: 2, Name: "humanNode3", Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	humanNode4 = &simpleGraph.Node{Id: 3, Name: "humanNode4", Seq: dna.StringToBases("CCC"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("CCC")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	humanNode5 = &simpleGraph.Node{Id: 4, Name: "humanNode5", Seq: dna.StringToBases("GGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("GGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}

	humanEdge1 = &simpleGraph.Edge{humanNode2, 0.25}
	humanEdge2 = &simpleGraph.Edge{humanNode3, 0.75}
	humanEdge3 = &simpleGraph.Edge{humanNode3, 1.00}
	humanEdge4 = &simpleGraph.Edge{humanNode4, 0.25}
	humanEdge5 = &simpleGraph.Edge{humanNode5, 0.75}
	humanEdge6 = &simpleGraph.Edge{humanNode5, 1.00}

	humanGraph = &simpleGraph.SimpleGraph{Nodes: []*simpleGraph.Node{humanNode1, humanNode2, humanNode3, humanNode4, humanNode5}}

	chimpNode1 = &simpleGraph.Node{Id: 0, Name: "chimpNode1", Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	chimpNode2 = &simpleGraph.Node{Id: 1, Name: "chimpNode2", Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}

	chimpEdge1 = &simpleGraph.Edge{chimpNode2, 1.00}

	chimpGraph = &simpleGraph.SimpleGraph{Nodes: []*simpleGraph.Node{chimpNode1, chimpNode2}}

	gorillaNode1 = &simpleGraph.Node{Id: 0, Name: "gorillaNode1", Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}
	gorillaNode2 = &simpleGraph.Node{Id: 1, Name: "gorillaNode2", Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil, Info: simpleGraph.Annotation{}}

	gorillaEdge1 = &simpleGraph.Edge{Dest: gorillaNode2, Prob: 1.00}

	gorillaGraph = &simpleGraph.SimpleGraph{Nodes: []*simpleGraph.Node{gorillaNode1, gorillaNode2}}

	nodeAlign0 = graphColumn{AlignId: 0, AlignNodes: map[string][]*simpleGraph.Node{"human": []*simpleGraph.Node{humanGraph.Nodes[0]}, "chimp": []*simpleGraph.Node{chimpGraph.Nodes[0]}, "gorilla": []*simpleGraph.Node{gorillaGraph.Nodes[0]}}}
	nodeAlign1 = graphColumn{AlignId: 1, AlignNodes: map[string][]*simpleGraph.Node{"human": []*simpleGraph.Node{humanGraph.Nodes[1]}}}
	nodeAlign2 = graphColumn{AlignId: 2, AlignNodes: map[string][]*simpleGraph.Node{"human": []*simpleGraph.Node{humanGraph.Nodes[2]}, "chimp": []*simpleGraph.Node{chimpGraph.Nodes[1]}, "gorilla": []*simpleGraph.Node{gorillaGraph.Nodes[1]}}}
	nodeAlign3 = graphColumn{AlignId: 3, AlignNodes: map[string][]*simpleGraph.Node{"human": []*simpleGraph.Node{humanGraph.Nodes[3], humanGraph.Nodes[4]}}}
)

func TestGraphColumn(t *testing.T) {
	allAlign = append(allAlign, nodeAlign0)
	allAlign = append(allAlign, nodeAlign1)
	allAlign = append(allAlign, nodeAlign2)
	allAlign = append(allAlign, nodeAlign3)

	for _, speciesNodes := range nodeAlign0.AlignNodes {
		if len(speciesNodes) != 1 {
			log.Fatal("incorrect number of human nodes")
		}
		_, ok := nodeAlign0.AlignNodes["human"]
		if !ok {
			log.Fatal("no human in align column")
		}
		_, ok = nodeAlign0.AlignNodes["chimp"]
		if !ok {
			log.Fatal("no chimp in align column")
		}
		_, ok = nodeAlign0.AlignNodes["gorilla"]
		if !ok {
			log.Fatal("no gorilla in align column")
		}
	}

	for _, speciesNodes := range nodeAlign1.AlignNodes {
		if len(speciesNodes) != 1 {
			log.Fatal("incorrect number of human nodes")
		}
		_, ok := nodeAlign1.AlignNodes["human"]
		if !ok {
			log.Fatal("no human in align column")
		}
		_, ok = nodeAlign1.AlignNodes["chimp"]
		if ok {
			log.Fatal("chimp in wrong align column")
		}
		_, ok = nodeAlign1.AlignNodes["gorilla"]
		if ok {
			log.Fatal("gorilla in wrong align column")
		}
	}

	for _, speciesNodes := range nodeAlign2.AlignNodes {
		if len(speciesNodes) != 1 {
			log.Fatal("incorrect number of human nodes")
		}
		_, ok := nodeAlign2.AlignNodes["human"]
		if !ok {
			log.Fatal("no human in align column")
		}
		_, ok = nodeAlign2.AlignNodes["chimp"]
		if !ok {
			log.Fatal("no chimp in align column")
		}
		_, ok = nodeAlign2.AlignNodes["gorilla"]
		if !ok {
			log.Fatal("no gorilla in align column")
		}
	}

	for _, speciesNodes := range nodeAlign3.AlignNodes {
		if len(speciesNodes) != 2 {
			log.Fatal("incorrect number of human nodes")
		}
		_, ok := nodeAlign3.AlignNodes["human"]
		if !ok {
			log.Fatal("no human in align column")
		}
		_, ok = nodeAlign3.AlignNodes["chimp"]
		if ok {
			log.Fatal("chimp in wrong align column")
		}
		_, ok = nodeAlign3.AlignNodes["gorilla"]
		if ok {
			log.Fatal("gorilla in wrong align column")
		}
	}

	humanNode1.Next = []*simpleGraph.Edge{humanEdge1, humanEdge2}
	humanNode2.Next = []*simpleGraph.Edge{humanEdge3}
	humanNode3.Next = []*simpleGraph.Edge{humanEdge4, humanEdge5}
	humanNode4.Next = []*simpleGraph.Edge{humanEdge6}

	humanNode2.Prev = []*simpleGraph.Edge{humanEdge1}
	humanNode3.Prev = []*simpleGraph.Edge{humanEdge2, humanEdge3}
	humanNode4.Prev = []*simpleGraph.Edge{humanEdge4}
	humanNode5.Prev = []*simpleGraph.Edge{humanEdge5, humanEdge6}

	chimpNode1.Next = []*simpleGraph.Edge{chimpEdge1}
	chimpNode2.Prev = []*simpleGraph.Edge{chimpEdge1}

	gorillaNode1.Next = []*simpleGraph.Edge{gorillaEdge1}
	gorillaNode2.Prev = []*simpleGraph.Edge{gorillaEdge1}

	//simpleGraph.PrintGraph(humanGraph)
	//simpleGraph.PrintGraph(chimpGraph)
	//simpleGraph.PrintGraph(gorillaGraph)
}

func TestPathFinder(t *testing.T) {
	actualSeq := "ACGTTTGGGGG"
	actualPath := []uint32{0, 2, 4}
	path, prob := PathFinder(humanGraph)
	for i := range path {
		if path[i] != actualPath[i] {
			log.Fatalf("path output: %v does not match expected path: %v", path, actualPath)
		}
	}
	if prob != 0.5625 {
		log.Fatal("prob of path is incorrect")
	}
	seq := seqOfPath(humanGraph, path)
	if dna.BasesToString(seq) != actualSeq {
		log.Printf("Sequence outpput: %s does not match expected sequence: %s", dna.BasesToString(seq), actualSeq)
	}
	//log.Print(dna.BasesToString(seq))
}

func TestBuildNodes(t *testing.T) {
	var id uint32
	var iNodes []*expandedTree.ETree
	tree, _ := expandedTree.ReadNewick("testdata/HCGAtree.txt")
	tNodes := expandedTree.GetTree(tree)

	for t := 0; t < len(tNodes); t++ {
		if tNodes[t].Right != nil && tNodes[t].Left != nil {
			iNodes = append(iNodes, tNodes[t])
		}
	}
	for in := 0; in < len(iNodes); in++ {
		speciesGraph := simpleGraph.NewGraph()
		for i := 0; i < len(allAlign); i++ {
			id = BuildNodes(iNodes[in], allAlign[i], id)
			for _, nodes := range allAlign[i].AlignNodes {
				for n := 0; n < len(nodes); n++ {
					if nodes[n].Name == iNodes[in].Name {
						speciesGraph.Nodes = append(speciesGraph.Nodes, nodes[n])
					}
				}
			}
		}
		if len(speciesGraph.Nodes) != 5 {
			log.Fatal("wrong number of nodes in the ancestor's graph")
		}
		//simpleGraph.PrintGraph(speciesGraph)
	}
}
