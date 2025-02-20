package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type Node struct {
	name        string
	ocr         bool
	score       int
	seen        bool
	connections []*Node
}

type bed struct {
	chrom   string
	start   string
	end     string
	names   []string
	famName []string
}

func main() {
	var cols, names, ocr, cp []string
	var tmpName string
	var j, k, homologous int
	var found bool
	var n *Node

	var startNode *string = flag.String("startNode", "", "Provide a node name to build a graph with. If \"all\" is used, a graph of all nodes will be created. The output will be directed to stdout.")
	flag.Parse()
	if len(flag.Args()) != 4 {
		fmt.Println("Usage: familiesFromLiftAndMerge in.merge.bed out.fam.bed out.fam.txt out.wholeMap.txt")
		log.Fatalf("Expected 4 args, got %d", len(flag.Args()))
	}

	mp := make(map[string]*Node)
	bedMap := make(map[string]bed)

	in := fileio.Read(flag.Arg(0))
	outBed := flag.Arg(1)
	outFamilies := flag.Arg(2)
	outMap := flag.Arg(3)

	for i := range in {
		cols = strings.Split(in[i], "\t")
		names = strings.Split(cols[3], ",")
		if len(names) == 1 && !strings.Contains(names[0], "lift") {
			bedMap[names[0]] = bed{chrom: cols[0], start: cols[1], end: cols[2], names: []string{names[0]}}
			continue
		}
		ocr = ocr[:0]
		for j = range names {
			if !strings.Contains(names[j], "lift") {
				ocr = append(ocr, names[j])
			}
		}
		if len(ocr) == 0 {
			ocr = append(ocr, fmt.Sprintf("homologous_%d", homologous))
			homologous++
		}

		for j = range ocr {
			cp = make([]string, len(ocr))
			copy(cp, ocr)
			bedMap[ocr[j]] = bed{chrom: cols[0], start: cols[1], end: cols[2], names: cp}
			_, found = mp[ocr[j]]
			if !found {
				createNode(mp, ocr[j])
			}
			for k = range names {
				if !strings.Contains(names[k], "lift") {
					continue
				}
				tmpName = stripLiftName(names[k])
				n, found = mp[tmpName]
				if !found {
					n = createNode(mp, names[k])
				}
				mp[ocr[j]].connections = checkConnectionSlice(mp[ocr[j]].connections, n)
				mp[n.name].connections = checkConnectionSlice(mp[n.name].connections, mp[ocr[j]])
			}
		}
	}

	//writeFamilies(mp, bedMap, outBed, flag.Arg(2))
	writeFamiliesRecursive(mp, bedMap, outFamilies, outBed)
	writeWholeMap(mp, outMap)
	if *startNode != "" {
		for i := range mp {
			mp[i].seen = false
		}
		buildGraphRecursive(mp, *startNode, "stdout")
	}
}

func createNode(mp map[string]*Node, enh string) *Node {
	var n *Node
	var openChrom bool = true

	if strings.Contains(enh, "homologous_") {
		openChrom = false
	}
	if strings.Contains(enh, "lift") {
		enh = stripLiftName(enh)
	}

	n = &Node{
		name:        enh,
		ocr:         openChrom,
		score:       0,
		seen:        false,
		connections: []*Node{},
	}
	mp[enh] = n
	return n
}

func writeFamiliesRecursive(mp map[string]*Node, bedMap map[string]bed, outFile string, outBedFile string) {
	var c int
	var fam []string
	out := fileio.EasyCreate(outFile)
	outBed := fileio.EasyCreate(outBedFile)
	for i := range mp {
		fam = fam[:0]
		if !mp[i].seen {
			c++
			fam = dfs(mp[i], fam)
		}
		if len(fam) > 0 {
			fileio.WriteToFileHandle(out, fmt.Sprintf("%d\tfamily_%d\t%s", len(fam), c, strings.Join(fam, ",")))
			addFamNameToBed(fam, bedMap, c)
		}
	}
	writeBedMap(bedMap, outBed)
	exception.PanicOnErr(out.Close())
}

func dfs(node *Node, fam []string) []string {
	if node.seen {
		return fam
	}

	node.seen = true
	fam = append(fam, node.name)

	for i := range node.connections {
		fam = dfs(node.connections[i], fam)
	}
	return fam
}

func checkConnectionSlice(conns []*Node, n *Node) []*Node {
	for i := range conns {
		if n.name == conns[i].name {
			return conns
		}
	}
	return append(conns, n)
}

func writeWholeMap(mp map[string]*Node, outfile string) {
	var names []string

	out := fileio.EasyCreate(outfile)

	for i := range mp {
		names = nodesToNameSlice(mp[i].connections)
		fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%s", i, strings.Join(names, ",")))
	}
	exception.PanicOnErr(out.Close())
}

//old code before recursive dfs implementation
/*func writeFamilies(mp map[string]*Node, bedMap map[string]bed, outBed *fileio.EasyWriter, outfile string) {
	var workSlice, currSlice []string
	var j, c int

	famMap := make(map[string]int)

	out := fileio.EasyCreate(outfile)

	for i := range mp {
		workSlice = []string{}
		currSlice = []string{}

		if mp[i].seen == true {
			continue
		}
		c++
		mp[i].seen = true

		famMap[i] = 0

		currSlice = nodesToNameSlice(mp[i].connections)
		addCurrSliceToMap(famMap, currSlice)

		workSlice = addToWorkSlice(workSlice, currSlice, mp)

		currSlice = []string{}

		for j = 0; j < len(workSlice); j++ {
			currSlice = nodesToNameSlice(mp[workSlice[j]].connections)
			addCurrSliceToMap(famMap, currSlice)
			workSlice = addToWorkSlice(workSlice, currSlice, mp)
			currSlice = []string{}
			mp[workSlice[j]].seen = true
		}

		addFamNameToBed(famMap, bedMap, c)
		writeFamMap(famMap, c, out)

		maps.Clear(famMap)

	}
	exception.PanicOnErr(out.Close())
	writeBedMap(bedMap, outBed)
}

*/

//for debug only
/*func printMap(mp map[string]*Node) {
	var conns []string
	for i := range mp {
		fmt.Println(mp[i])
		conns = nodesToNameSlice(mp[i].connections)
		for j := range conns {
			mp[i].connections[j].score = 7
			fmt.Printf("\t%v\n\t%v\n", mp[conns[j]], mp[i].connections[j])
		}
	}
}

*/

func nodesToNameSlice(nodeSlice []*Node) []string {
	var stringSlice []string

	for i := range nodeSlice {
		stringSlice = append(stringSlice, nodeSlice[i].name)
	}
	return stringSlice
}

//part of old family partition code
/*func addCurrSliceToMap(famMap map[string]int, currSlice []string) {
	for i := range currSlice {
		famMap[currSlice[i]] = 0
	}
}

*/

func writeBedMap(bedMap map[string]bed, outBed *fileio.EasyWriter) {
	var b bed
	for i := range bedMap {
		b, _ = bedMap[i]
		if len(b.famName) == 0 {
			b.famName = []string{"lonely"}
		}
		fileio.WriteToFileHandle(outBed, fmt.Sprintf("%s\t%s\t%s\t%s\t%s", b.chrom, b.start, b.end, strings.Join(b.names, ","), strings.Join(b.famName, ",")))
	}
	exception.PanicOnErr(outBed.Close())
}

func addFamNameToBed(famSlice []string, bedMap map[string]bed, c int) {
	var found bool
	var b bed
	for i := range famSlice {
		b, found = bedMap[famSlice[i]]
		if !found {
			fmt.Println("bed record not found xyz ", i)
			continue
		}
		b.famName = append(b.famName, fmt.Sprintf("family_%d", c))
		bedMap[famSlice[i]] = b
	}
}

// part of old family partition code
func addToWorkSlice(workSlice, currSlice []string, mp map[string]*Node) []string {
	for i := range currSlice {
		if !mp[currSlice[i]].seen {
			workSlice = append(workSlice, currSlice[i])
		}
	}
	return workSlice
}

//part of old family partition code
/*func writeFamMap(famMap map[string]int, c int, out *fileio.EasyWriter) {
	var names []string
	for i := range famMap {
		names = append(names, i)
	}
	fileio.WriteToFileHandle(out, fmt.Sprintf("%d\tfamily_%d\t%s", len(names), c, strings.Join(names, ",")))
}

*/

func stripLiftName(name string) string {
	return strings.Join(strings.Split(name, "_")[:5], "_")
}

func buildGraphRecursive(mp map[string]*Node, startNode string, outDot string) {
	out := fileio.EasyCreate(outDot)
	fileio.WriteToFileHandle(out, "strict graph {")

	if startNode == "all" {
		for i := range mp {
			dfsDotFile(mp[i], "", out)
		}
		fileio.WriteToFileHandle(out, "}")
		exception.PanicOnErr(out.Close())
		return
	}

	fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [fillcolor=green, style=filled, label=\"\"]", startNode))
	n := mp[startNode]

	dfsDotFile(n, startNode, out)

	fileio.WriteToFileHandle(out, "}")
	exception.PanicOnErr(out.Close())
}

func dfsDotFile(node *Node, startNode string, out *fileio.EasyWriter) {
	if node.seen {
		return
	}
	node.seen = true
	if node.ocr && node.name != startNode {
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [fillcolor=red, style=filled, label=\"\"]", node.name))
	}
	connections, homologous := formatAlnString(node.connections)
	for i := range homologous {
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [label=\"\"]", homologous[i]))
	}
	fileio.WriteToFileHandle(out, fmt.Sprintf("\t\"%s\" -- %s", node.name, connections))

	for i := range node.connections {
		dfsDotFile(node.connections[i], startNode, out)
	}
}

/*func buildGraph(mp map[string]*Node, startNode string, outDot string) {
	var workSlice, currSlice []string

	if startNode == "all" {
		//buildGraphAllNodes(mp, outDot)
	} else {
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

}

*/

func formatAlnString(nodeSlice []*Node) (string, []string) {
	var sb strings.Builder
	var homologous []string
	sb.WriteString("{")
	for i := range nodeSlice {
		sb.WriteString("\"" + nodeSlice[i].name + "\"" + " ")
		if strings.Contains(nodeSlice[i].name, "homologous_") {
			homologous = append(homologous, nodeSlice[i].name)
		}
	}
	sb.WriteString("}")
	return sb.String(), homologous
}
