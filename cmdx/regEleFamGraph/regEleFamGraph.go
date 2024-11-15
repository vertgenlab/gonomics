package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"strings"
)

type Node struct {
	name        string
	coord       bed.Bed
	ocr         bool
	score       int
	seen        bool
	connections []*Node
}

func main() {
	var cols, alns []string
	var n *Node
	var found, ocr bool
	var name, alnName string
	var b bed.Bed
	var conns []*Node

	flag.Parse()

	if len(flag.Args()) != 3 {
		log.Fatalf("~/go/bin/regEleFamGraph in.txt centralNode out.dot")
	}

	mp := make(map[string]*Node)
	txt := fileio.Read(flag.Arg(0))

	for i := 0; i < len(txt); i++ {
		cols = strings.Split(txt[i], "\t")
		if cols[5] == "" {
			continue
		}
		name = cols[3]
		n, found = mp[cols[3]]
		if !found {
			n = &Node{
				name:        name,
				coord:       bed.Bed{Chrom: cols[0], ChromStart: parse.StringToInt(cols[1]), ChromEnd: parse.StringToInt(cols[2]), FieldsInitialized: 3},
				ocr:         true,
				seen:        false,
				connections: []*Node{},
				score:       parse.StringToInt(cols[4]),
			}
			mp[name] = n
		} else {
			n.score = parse.StringToInt(cols[4])
		}

		conns = mp[name].connections

		alns = parseAlns(cols[5])
		for j := range alns {
			alnName, b, ocr = parseSingleAln(alns[j])
			n, found = mp[alnName]
			if found {
				conns = append(conns, n)
			} else {
				n = &Node{
					name:        alnName,
					coord:       b,
					ocr:         ocr,
					seen:        false,
					connections: []*Node{},
				}
				conns = append(conns, n)
				mp[alnName] = n
			}
			if strings.Contains(n.name, "lift") {
				n.connections = append(n.connections, mp[name])
			}
		}
		mp[name].connections = conns
	}
	buildGraph(mp, flag.Arg(1), flag.Arg(2))
}

func parseAlns(s string) []string {
	slc := strings.Split(s, ";")
	return slc[:len(slc)-1]
}

func parseSingleAln(s string) (name string, b bed.Bed, ocr bool) {
	slc := strings.Split(s, "_")
	b = bed.Bed{Chrom: slc[0], ChromStart: parse.StringToInt(slc[1]), ChromEnd: parse.StringToInt(slc[2]), FieldsInitialized: 3}
	name = strings.Join(slc[3:], "_")
	if strings.Contains(name, "lift") {
		ocr = false
	} else {
		ocr = true
	}
	return name, b, ocr
}

func buildGraph(mp map[string]*Node, startNode string, outDot string) {
	var workSlice, currSlice []string

	if startNode == "all" {
		buildGraphAllNodes(mp, outDot)
	}

	out := fileio.EasyCreate(outDot)
	fileio.WriteToFileHandle(out, "strict graph {")

	currSlice = nodesToNameSlice(mp[startNode].connections)
	fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [fillcolor=green, style=filled]", startNode))
	fileio.WriteToFileHandle(out, fmt.Sprintf("\t\"%s\" -- %s", startNode, formatAlnString(currSlice)))
	workSlice = addToWorkSlice(workSlice, currSlice, mp)
	mp[startNode].seen = true

	for i := 0; i < len(workSlice); i++ {
		//fmt.Println(mp[workSlice[i]])
		currSlice = nodesToNameSlice(mp[workSlice[i]].connections)
		//fmt.Println(currSlice)
		workSlice = addToWorkSlice(workSlice, currSlice, mp)
		if !strings.Contains(workSlice[i], "lift") {
			fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [fillcolor=red, style=filled]", workSlice[i]))
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t\"%s\" -- %s", workSlice[i], formatAlnString(currSlice)))
		mp[workSlice[i]].seen = true
	}

	fileio.WriteToFileHandle(out, "}")
	exception.PanicOnErr(out.Close())

}

func formatAlnString(currSlice []string) string {
	var sb strings.Builder
	sb.WriteString("{")
	for i := range currSlice {
		sb.WriteString("\"" + currSlice[i] + "\"" + " ")
	}
	sb.WriteString("}")
	return sb.String()
}

func nodesToNameSlice(nodeSlice []*Node) []string {
	var stringSlice []string

	for i := range nodeSlice {
		stringSlice = append(stringSlice, nodeSlice[i].name)
	}
	return stringSlice
}

func addToWorkSlice(workSlice, currSlice []string, mp map[string]*Node) []string {
	for i := range currSlice {
		if !mp[currSlice[i]].seen {
			workSlice = append(workSlice, currSlice[i])
		}
	}
	return workSlice
}

func buildGraphAllNodes(mp map[string]*Node, outDot string) {

	out := fileio.EasyCreate(outDot)
	fileio.WriteToFileHandle(out, "strict graph {")

	for i := range mp {
		if !strings.Contains(i, "lift") {
			fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [fillcolor=red, style=filled]", i))
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t\"%s\" -- %s", i, formatAlnString(nodesToNameSlice(mp[i].connections))))
	}
	fileio.WriteToFileHandle(out, "}")
	exception.PanicOnErr(out.Close())
}
