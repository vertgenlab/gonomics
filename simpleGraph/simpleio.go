package simpleGraph

import (
	"bytes"
	"github.com/edotau/simpleio"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"strings"
	"sync"
)

func NewGenomeGraph() SimpleGraph {
	graph := new(SimpleGraph)
	graph.Nodes = make([]*Node, 0)
	return *graph
}

func SimplyRead(filename string) SimpleGraph {
	simpleioReader := simpleio.NewSimpleReader(filename)
	genome := NewGenomeGraph()
	var currNode *Node
	var edges map[string]*Node = make(map[string]*Node)
	var weight float32

	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	var data *bytes.Buffer
	var line string
	var words []string = make([]string, 0, 2)
	var text []string = make([]string, 0, 3)
	var seqIdx int32 = -1
	var i int
	var ok bool

	for reader, done := simpleio.ReadLine(simpleioReader); !done; reader, done = simpleio.ReadLine(simpleioReader) {
		data = simplePool.Get().(*bytes.Buffer)
		data.Write(reader)
		line = data.String()
		switch true {
		case strings.HasPrefix(line, ">"):
			seqIdx++
			words = strings.Split(line[1:], ":")

			currNode = &Node{Id: uint32(seqIdx), Name: string(words[0])}
			if len(words) == 2 {
				text = strings.Split(words[1], "_")
				currNode.Info = &Annotation{Allele: uint8(common.StringToUint32(text[1])), Start: common.StringToUint32(text[3]), Variant: uint8(common.StringToUint32(text[2]))}
			}
			AddNode(&genome, currNode)
			_, ok = edges[line[1:]]
			if !ok {
				edges[string(line[1:])] = currNode
			}
		case strings.Contains(line, "\t"):
			words = strings.Split(line, "\t")
			if len(words) > 2 {
				for i = 1; i < len(words); i += 2 {
					weight = float32(common.StringToFloat64(words[i]))
					AddEdge(edges[words[0]], edges[words[i+1]], weight)
				}
			}
		case !strings.ContainsAny(line, "\t:"):
			genome.Nodes[seqIdx].Seq = append(genome.Nodes[seqIdx].Seq, simpleio.ByteSliceToDnaBases(data.Bytes())...)
		}
		data.Reset()
		simplePool.Put(data)
	}
	simpleioReader.Close()
	for i = 0; i < len(genome.Nodes); i++ {
		genome.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(genome.Nodes[i].Seq)
	}
	return genome
}

func readFastqGsw(fileOne string, fileTwo string, answer chan<- *fastq.PairedEndBig) {
	readOne, readTwo := simpleio.NewSimpleReader(fileOne), simpleio.NewSimpleReader(fileTwo)
	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	for curr, done := fqPair(readOne, readTwo, simplePool); !done; curr, done = fqPair(readOne, readTwo, simplePool) {
		answer <- curr
	}
	close(answer)
}

func fqPair(reader1 *simpleio.SimpleReader, reader2 *simpleio.SimpleReader, simplePool sync.Pool) (*fastq.PairedEndBig, bool) {
	fqOne, done1 := nextFq(reader1, simplePool)
	fqTwo, done2 := nextFq(reader2, simplePool)
	if (!done1 && done2) || (done1 && !done2) {
		log.Fatalf("Error: fastq files do not end at the same time...\n")
	} else if done1 || done2 {
		return nil, true
	}
	curr := NewBigFastqPair()
	curr.Fwd, curr.Rev = fqOne, fqTwo

	//curr.Fwd.Name, curr.Rev.Name = strings.Split(fqOne.Name, " ")[0], strings.Split(fqTwo.Name, " ")[0]
	return curr, false
}

func nextFq(reader *simpleio.SimpleReader, simplePool sync.Pool) (*fastq.FastqBig, bool) {
	name, done := simpleio.ReadLine(reader)
	if done {
		return nil, true
	}
	seq, sDone := simpleio.ReadLine(reader)
	plus, pDone := simpleio.ReadLine(reader)
	qual, qDone := simpleio.ReadLine(reader)

	if sDone || pDone || qDone {
		log.Fatalf("Error: There is an empty line in this fastq record\n")
	}

	answer := fastq.FastqBig{}
	data := simplePool.Get().(*bytes.Buffer)
	data.Write(name)
	answer.Name = strings.Split(data.String()[1:], " ")[0]
	data.Reset()
	data.Write(seq)
	//set up sequence and reverse comp
	answer.Seq = simpleio.ByteSliceToDnaBases(data.Bytes())
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)

	//performs two bit conversion
	answer.Rainbow = dnaTwoBit.NewTwoBitRainbow(answer.Seq)
	answer.RainbowRc = dnaTwoBit.NewTwoBitRainbow(answer.SeqRc)

	data.Reset()

	if string(plus) != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}

	data.Write(qual)
	answer.Qual = fastq.ToQualUint8(bytes.Runes(data.Bytes()))
	data.Reset()
	simplePool.Put(data)
	return &answer, false
}

func NewBigFastqPair() *fastq.PairedEndBig {
	return &fastq.PairedEndBig{
		Fwd: new(fastq.FastqBig),
		Rev: new(fastq.FastqBig),
	}
}
