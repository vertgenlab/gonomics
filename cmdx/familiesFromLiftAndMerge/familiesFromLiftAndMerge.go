package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"golang.org/x/exp/maps"
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
	var j, k, homologous, multi int
	var found bool
	var n *Node

	flag.Parse()

	mp := make(map[string]*Node)
	bedMap := make(map[string]bed)

	in := fileio.Read(flag.Arg(0))
	outBed := fileio.EasyCreate(flag.Arg(1))

	for i := range in {
		cols = strings.Split(in[i], "\t")
		names = strings.Split(cols[3], ",")
		if len(names) == 1 && !strings.Contains(names[0], "lift") {
			//fileio.WriteToFileHandle(outBed, strings.Join(cols[0:4], "\t"))
			continue
		}
		ocr = ocr[:0]
		for j = range names {
			if !strings.Contains(names[j], "lift") {
				ocr = append(ocr, names[j])
			}
		}
		switch len(ocr) {
		case 0:
			ocr = append(ocr, fmt.Sprintf("homologous_%d", homologous))
			homologous++
		case 1:
		default:
			fmt.Println(ocr)
			multi++
		}
		for j = range ocr {
			cp = make([]string, len(ocr))
			copy(cp, ocr)
			bedMap[ocr[j]] = bed{chrom: cols[0], start: cols[1], end: cols[2], names: cp}
			createNode(mp, ocr[j])
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
	fmt.Println(multi)
	writeFamilies(mp, bedMap, outBed, flag.Arg(2))
	//writeWholeMap(mp, flag.Arg(3))
	//buildGraph(mp, "h9_atac_rep1_peak_14853", "stdout")
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
	var conns []*Node
	var names []string
	var j int

	out := fileio.EasyCreate(outfile)

	for i := range mp {
		conns = conns[:0]
		conns = mp[i].connections
		for j = range conns {
			names = append(names, conns[j].name)
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%s", i, strings.Join(names, ",")))
	}
	exception.PanicOnErr(out.Close())

}

func writeFamilies(mp map[string]*Node, bedMap map[string]bed, outBed *fileio.EasyWriter, outfile string) {
	var workSlice, currSlice []string
	var j, c int

	famMap := make(map[string]int)

	out := fileio.EasyCreate(outfile)

	for i := range mp {
		workSlice = workSlice[:0]
		currSlice = currSlice[:0]

		if mp[i].seen == true {
			continue
		}
		c++
		mp[i].seen = true

		famMap[i] = 0

		currSlice = nodesToNameSlice(mp[i].connections)
		addCurrSliceToMap(famMap, currSlice)

		workSlice = addToWorkSlice(workSlice, currSlice, mp)

		for j = 0; j < len(workSlice); j++ {
			currSlice = nodesToNameSlice(mp[workSlice[j]].connections)
			addCurrSliceToMap(famMap, currSlice)
			workSlice = addToWorkSlice(workSlice, currSlice, mp)
			mp[workSlice[j]].seen = true
		}

		addFamNameToBed(famMap, bedMap, c)
		writeFamMap(famMap, c, out)

		maps.Clear(famMap)

	}
	exception.PanicOnErr(out.Close())
	writeBedMap(bedMap, outBed)
}

func writeBedMap(bedMap map[string]bed, outBed *fileio.EasyWriter) {
	var b bed
	for i := range bedMap {
		b, _ = bedMap[i]
		fileio.WriteToFileHandle(outBed, fmt.Sprintf("%s\t%s\t%s\t%s\t%s", b.chrom, b.start, b.end, strings.Join(b.names, ","), strings.Join(b.famName, ",")))
	}
	exception.PanicOnErr(outBed.Close())
}

func addFamNameToBed(famMap map[string]int, bedMap map[string]bed, c int) {
	var found bool
	var b bed
	for i := range famMap {
		b, found = bedMap[i]
		if !found {
			continue
		}
		b.famName = append(b.famName, fmt.Sprintf("family_%d", c))
		bedMap[i] = b
	}
}

func nodesToNameSlice(nodeSlice []*Node) []string {
	var stringSlice []string

	for i := range nodeSlice {
		stringSlice = append(stringSlice, nodeSlice[i].name)
	}
	return stringSlice
}

func addCurrSliceToMap(famMap map[string]int, currSlice []string) {
	for i := range currSlice {
		famMap[currSlice[i]] = 0
	}
}

func addToWorkSlice(workSlice, currSlice []string, mp map[string]*Node) []string {
	for i := range currSlice {
		if !mp[currSlice[i]].seen {
			workSlice = append(workSlice, currSlice[i])
		}
	}
	return workSlice
}

func writeFamMap(famMap map[string]int, c int, out *fileio.EasyWriter) {
	var names []string
	for i := range famMap {
		names = append(names, i)
	}
	fileio.WriteToFileHandle(out, fmt.Sprintf("%d\tfamily_%d\t%s", len(names), c, strings.Join(names, ",")))
}

func createNode(mp map[string]*Node, enh string) *Node {
	var n *Node
	var openChrom bool = true

	if strings.Contains(enh, "homologous_") {
		openChrom = false
	}
	if strings.Contains(enh, "lift") {
		openChrom = false
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

func stripLiftName(name string) string {
	return strings.Join(strings.Split(name, "_")[:5], "_")
}

func buildGraph(mp map[string]*Node, startNode string, outDot string) {
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

func formatAlnString(currSlice []string) string {
	var sb strings.Builder
	sb.WriteString("{")
	for i := range currSlice {
		sb.WriteString("\"" + currSlice[i] + "\"" + " ")
	}
	sb.WriteString("}")
	return sb.String()
}
