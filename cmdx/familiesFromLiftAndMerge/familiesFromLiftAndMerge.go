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
	percID      []float64
}

type bed struct {
	chrom     string
	start     string
	end       string
	names     []string
	famName   []string
	PercID    []float64
	writeCopy bool
}

func main() {
	var cols, names, ocr, cp []string
	var percID []float64
	var tmpName string
	var j, k, d int
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
		percID = fileio.StringToFloatSlice(cols[6])
		if len(names) == 1 && !strings.Contains(names[0], "lift") {
			bedMap[names[0]] = bed{chrom: cols[0], start: cols[1], end: cols[2], names: []string{names[0]}, PercID: percID, writeCopy: true}
			continue
		}
		ocr = ocr[:0]
		d = 0

		for j = range names {
			if !strings.Contains(names[j], "lift") {
				ocr = append(ocr, names[j])
				d++
			}
		}
		for j = range ocr {
			cp = make([]string, len(ocr))
			copy(cp, ocr)
			if j == 0 {
				bedMap[ocr[j]] = bed{chrom: cols[0], start: cols[1], end: cols[2], PercID: percID, names: cp, writeCopy: true}
			} else {
				bedMap[ocr[j]] = bed{chrom: cols[0], start: cols[1], end: cols[2], PercID: percID, names: cp, writeCopy: false}
			}

			_, found = mp[ocr[j]]
			if !found {
				createNode(mp, ocr[j])
			}
			for k = range names {
				if !strings.Contains(names[k], "lift") || strings.Contains(names[k], ocr[j]) {
					continue
				}

				tmpName = stripLiftName(names[k])
				n, found = mp[tmpName]
				if !found {
					n = createNode(mp, names[k])
				}
				if addConnectionSlice(mp[ocr[j]].connections, n) {
					mp[ocr[j]].connections = append(mp[ocr[j]].connections, n)
					mp[n.name].percID = append(mp[n.name].percID, percID[k])
				} else {
					mp[n.name].percID = append(mp[n.name].percID, percID[k])
				}
				if addConnectionSlice(mp[n.name].connections, mp[ocr[j]]) {
					mp[n.name].connections = append(mp[n.name].connections, mp[ocr[j]])
				}
			}
		}
	}

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

	if strings.Contains(enh, "homologous") {
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
		percID:      []float64{},
	}
	mp[enh] = n
	return n
}

func writeFamiliesRecursive(mp map[string]*Node, bedMap map[string]bed, outFile string, outBedFile string) {
	var c int
	var fam []string
	var pidSlice []float64
	var min, max, avg float64
	out := fileio.EasyCreate(outFile)
	outBed := fileio.EasyCreate(outBedFile)
	for i := range mp {
		fam = fam[:0]
		pidSlice = pidSlice[:0]
		if !mp[i].seen {
			pidSlice = append(pidSlice, mp[i].percID[0])
			c++
			fam = dfs(mp[i], fam, pidSlice)
		}
		if len(fam) > 0 {
			min, max, avg = percIdCalcs(pidSlice)
			fileio.WriteToFileHandle(out, fmt.Sprintf("%d\tfamily_%d\t%s\t%.2f\t%.2f\t%.2f", len(fam), c, strings.Join(fam, ","), avg, min, max))
			addFamNameToBed(fam, bedMap, c)
		}
	}
	writeBedMap(bedMap, outBed)
	exception.PanicOnErr(out.Close())
}

func dfs(node *Node, fam []string, pidSlice []float64) []string {
	if node.seen {
		return fam
	}

	node.seen = true
	fam = append(fam, node.name)

	for i := range node.connections {
		pidSlice = append(pidSlice, node.percID[i])
		fam = dfs(node.connections[i], fam, pidSlice)
	}
	return fam
}

func percIdCalcs(pidSlice []float64) (min, max, avg float64) {
	var tot float64
	min, max = pidSlice[0], pidSlice[0]
	for i := range pidSlice {
		if pidSlice[i] > max {
			max = pidSlice[i]
		}
		if pidSlice[i] < min {
			min = pidSlice[i]
		}
		tot += pidSlice[i]
	}
	return min, max, tot / float64(len(pidSlice))
}

func addConnectionSlice(conns []*Node, n *Node) bool {
	for i := range conns {
		if n.name == conns[i].name {
			return false
		}
	}
	return true
}

func writeWholeMap(mp map[string]*Node, outfile string) {
	var names []string

	out := fileio.EasyCreate(outfile)

	for i := range mp {
		names = nodesToNameSlice(mp[i].connections)
		fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%s\t%s", i, strings.Join(names, ","), fileio.FloatSliceToString(mp[i].percID)))
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
		if !b.writeCopy {
			continue
		}
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
/*
func addToWorkSlice(workSlice, currSlice []string, mp map[string]*Node) []string {
	for i := range currSlice {
		if !mp[currSlice[i]].seen {
			workSlice = append(workSlice, currSlice[i])
		}
	}
	return workSlice
}

//part of old family partition code
func writeFamMap(famMap map[string]int, c int, out *fileio.EasyWriter) {
	var names []string
	for i := range famMap {
		names = append(names, i)
	}
	fileio.WriteToFileHandle(out, fmt.Sprintf("%d\tfamily_%d\t%s", len(names), c, strings.Join(names, ",")))
}

*/

func stripLiftName(name string) string {
	tmp := strings.Split(name, "_")
	return strings.Join(tmp[:len(tmp)-1], "_")
}

func buildGraphRecursive(mp map[string]*Node, startNode string, outDot string) {
	out := fileio.EasyCreate(outDot)
	fileio.WriteToFileHandle(out, "strict graph {")
	fileio.WriteToFileHandle(out, fmt.Sprintf("outputorder=\"edgesfirst\""))
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

	for i := range node.connections {
		if strings.Contains(node.connections[i].name, "homologous") && node.connections[i].name != startNode {
			fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [fillcolor=white, style=filled, label=\"\"]", node.connections[i].name))
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t\"%s\" -- \"%s\" [minlen = %.2f];", node.name, node.connections[i].name, 100-node.percID[i]))
	}

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
		if strings.Contains(nodeSlice[i].name, "homologous") {
			homologous = append(homologous, nodeSlice[i].name)
		}
	}
	sb.WriteString("}")
	return sb.String(), homologous
}
