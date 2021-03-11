package genomeGraph

import (
//"fmt"
//"github.com/vertgenlab/gonomics/common"
//"github.com/vertgenlab/gonomics/dna"
//"github.com/vertgenlab/gonomics/fileio"
//"log"
//"strings"
)

//TODO: work in progress
/*func DotToGraph(input string) *SimpleGraph {
	var line string
	var gg *SimpleGraph = NewGraph()
	var doneReading bool = false

	file := fileio.EasyOpen(input)
	defer file.Close()
	var nodeId uint32 = 0
	var words []string

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.Contains(line, "digraph") {
			log.Printf("Header line: digraph\n")
		} else if !strings.Contains(line, " -> ") && !strings.Contains(line, "//") && !strings.Contains(line, "}") {
			words = strings.Split(line, " ")
			AddNode(gg, &Node{Id: nodeId, Name: fmt.Sprint(nodeId)})

			if strings.Contains(line, "label") {
				workSpace := strings.Split(words[1], "=")
				//log.Printf("%v\n", strings.Trim(workSpace[1], "]"))
				gg.Nodes[nodeId].Seq = dna.StringToBases(strings.Trim(workSpace[1], "]"))
			}
			nodeId++
		} else if strings.Contains(line, " -> ") && !strings.Contains(line, "//") {
			words = strings.Split(line, "\t")
			//log.Printf("%v\n", words[1])
			workSpace := strings.Split(words[1], " -> ")
			AddEdge(gg.Nodes[common.StringToUint32(workSpace[0])], gg.Nodes[common.StringToUint32(workSpace[1])], 1)
		}
	}
	return gg
}*/
