package graphReconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"log"
	"testing"
)

var (
	allAlign []graphColumn

	humanGraph = genomeGraph.EmptyGraph()

	humanNode1 = genomeGraph.AddNode(humanGraph, &genomeGraph.Node{Id: 0, Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil})
	humanNode2 = genomeGraph.AddNode(humanGraph, &genomeGraph.Node{Id: 1, Seq: dna.StringToBases("AAA"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("AAA")), Prev: nil, Next: nil})
	humanNode3 = genomeGraph.AddNode(humanGraph, &genomeGraph.Node{Id: 2, Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil})
	humanNode4 = genomeGraph.AddNode(humanGraph, &genomeGraph.Node{Id: 3, Seq: dna.StringToBases("CCC"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("CCC")), Prev: nil, Next: nil})
	humanNode5 = genomeGraph.AddNode(humanGraph, &genomeGraph.Node{Id: 4, Seq: dna.StringToBases("GGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("GGG")), Prev: nil, Next: nil})

	humanEdge1 = genomeGraph.Edge{humanNode2, 0.25}
	humanEdge2 = genomeGraph.Edge{humanNode3, 0.75}
	humanEdge3 = genomeGraph.Edge{humanNode3, 1.00}
	humanEdge4 = genomeGraph.Edge{humanNode4, 0.25}
	humanEdge5 = genomeGraph.Edge{humanNode5, 0.75}
	humanEdge6 = genomeGraph.Edge{humanNode5, 1.00}

	chimpGraph = genomeGraph.EmptyGraph()

	chimpNode1 = genomeGraph.AddNode(chimpGraph, &genomeGraph.Node{Id: 0, Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil})
	chimpNode2 = genomeGraph.AddNode(chimpGraph, &genomeGraph.Node{Id: 1, Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil})

	chimpEdge1 = genomeGraph.Edge{chimpNode2, 1.00}

	gorillaGraph = genomeGraph.EmptyGraph()

	gorillaNode1 = genomeGraph.AddNode(gorillaGraph, &genomeGraph.Node{Id: 0, Seq: dna.StringToBases("ACGT"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ACGT")), Prev: nil, Next: nil})
	gorillaNode2 = genomeGraph.AddNode(gorillaGraph, &genomeGraph.Node{Id: 1, Seq: dna.StringToBases("TTGG"), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("TTGG")), Prev: nil, Next: nil})

	gorillaEdge1 = genomeGraph.Edge{Dest: gorillaNode2, Prob: 1.00}

	nodeAlign0 = graphColumn{AlignId: 0, AlignNodes: map[string][]*genomeGraph.Node{"human": []*genomeGraph.Node{&humanGraph.Nodes[0]}, "chimp": []*genomeGraph.Node{&chimpGraph.Nodes[0]}, "gorilla": []*genomeGraph.Node{&gorillaGraph.Nodes[0]}}}
	nodeAlign1 = graphColumn{AlignId: 1, AlignNodes: map[string][]*genomeGraph.Node{"human": []*genomeGraph.Node{&humanGraph.Nodes[1]}}}
	nodeAlign2 = graphColumn{AlignId: 2, AlignNodes: map[string][]*genomeGraph.Node{"human": []*genomeGraph.Node{&humanGraph.Nodes[2]}, "chimp": []*genomeGraph.Node{&chimpGraph.Nodes[1]}, "gorilla": []*genomeGraph.Node{&gorillaGraph.Nodes[1]}}}
	nodeAlign3 = graphColumn{AlignId: 3, AlignNodes: map[string][]*genomeGraph.Node{"human": []*genomeGraph.Node{&humanGraph.Nodes[3], &humanGraph.Nodes[4]}}}
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

	humanNode1.Next = []genomeGraph.Edge{humanEdge1, humanEdge2}
	humanNode2.Next = []genomeGraph.Edge{humanEdge3}
	humanNode3.Next = []genomeGraph.Edge{humanEdge4, humanEdge5}
	humanNode4.Next = []genomeGraph.Edge{humanEdge6}

	humanNode2.Prev = []genomeGraph.Edge{humanEdge1}
	humanNode3.Prev = []genomeGraph.Edge{humanEdge2, humanEdge3}
	humanNode4.Prev = []genomeGraph.Edge{humanEdge4}
	humanNode5.Prev = []genomeGraph.Edge{humanEdge5, humanEdge6}

	chimpNode1.Next = []genomeGraph.Edge{chimpEdge1}
	chimpNode2.Prev = []genomeGraph.Edge{chimpEdge1}

	gorillaNode1.Next = []genomeGraph.Edge{gorillaEdge1}
	gorillaNode2.Prev = []genomeGraph.Edge{gorillaEdge1}

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
	var treeFilename string = "testdata/HCGAtree.newick"
	var id uint32
	var speciesNodeCount int
	var internalNodes []*expandedTree.ETree
	tree, err := expandedTree.ReadNewick(treeFilename)
	if err != nil {
		t.Errorf("Error: unable to read in %s.  Got error: %s\n", treeFilename, err)
	}
	treeNodes := expandedTree.GetTree(tree)

	for t := 0; t < len(treeNodes); t++ {
		if treeNodes[t].Right != nil && treeNodes[t].Left != nil {
			internalNodes = append(internalNodes, treeNodes[t])
		}
	}
	for in := 0; in < len(internalNodes); in++ {
		speciesNodeCount = 0
		for i := 0; i < len(allAlign); i++ {
			id = BuildNodes(internalNodes[in], allAlign[i], id)
			nodesForAncestorAtCol := allAlign[i].AlignNodes[internalNodes[in].Name]
			speciesNodeCount += len(nodesForAncestorAtCol)
		}
		if speciesNodeCount != 5 {
			t.Error("Error: wrong number of nodes in the ancestor's graph\n")
		}
	}
}
