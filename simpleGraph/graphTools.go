package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"strings"
	//"runtime"
	//"sync"
	"os"
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
			
			if vcfFlag == 0 {
				currMatch.Seq = currMatch.Seq[:len(currMatch.Seq)-1]
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

func GSWAlignerOne(ref *SimpleGraph, input string, output string) {
	var tileSize int = 12
	var stepSize int = 1

	//wg := new(sync.WaitGroup)
	outFile, _ := os.Create(output)

	defer outFile.Close()

	header := NodesHeader(ref.Nodes)


	sam.WriteHeaderToFileHandle(outFile, header)
	file := fileio.EasyOpen(input)
	defer file.Close()
	//runtime.GOMAXPROCS(runtime.NumCPU())

	//threads := runtime.NumCPU()
	//tasks := make(chan *fastq.Fastq, 20)

	log.Printf("Indexing the genome...\n")
	ham5 := IndexGenomeDev(ref.Nodes, tileSize, stepSize)
	m, trace := SwMatrixSetup(10000)
	log.Printf("Aligning reads...\n")

	//var mappedRead *sam.SamAln
	/*
	for workers := 0; workers < threads; workers++ {
		go func(jobNumber int) {
			defer wg.Done()
			for {
				query, ok := <-tasks
				if !ok {
					return
				}
				mappedRead = GraphSmithWaterman(ref, query, ham5, tileSize, m, trace)
				sam.WriteAlnToFileHandle(outFile, mappedRead)
				//MapFastq(ref, query, seedLen, m, mSet, trace, outFile)
			}

		}(workers)
	}*/

	var fq *fastq.Fastq
	var done bool
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		//mappedRead = GraphSmithWaterman(ref, fq, ham5, tileSize, m, trace)
		//sam.WriteAlnToFileHandle(outFile, mappedRead)
		//wg.Add(1)
		//tasks <- fq
		go routinesGenomeGraph(ref, fq, ham5, tileSize, m, trace, outFile)
	}
	//close(tasks)
	//wg.Wait()
}
//Function calls GSW alignment, meant to be used in goroutines, and writes the alignment stgraight to sam file
func routinesGenomeGraph(gg *SimpleGraph, read *fastq.Fastq, seedHash [][]*SeedBed, seedLen int, m [][]int64, trace [][]rune, output *os.File) {
	var mappedRead *sam.SamAln
	mappedRead = GraphSmithWaterman(gg, read, seedHash, seedLen, m, trace)
	sam.WriteAlnToFileHandle(output, mappedRead)
}
