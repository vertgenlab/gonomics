package simpleGraph

import(
	"github.com/edotau/simpleio"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"bytes"
	"sync"
	"strings"

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



func ReadFastqGsw() {


}