package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

func main() {
	var c int
	var cols []string
	var name, connections string
	flag.Parse()
	out := fileio.EasyCreate(flag.Arg(1))
	txt := fileio.Read(flag.Arg(0))
	fileio.WriteToFileHandle(out, "graph regEleFam {")
	for i := range txt {
		cols = strings.Split(txt[i], "\t")
		if cols[5] == "" {
			continue
		}
		c++
		name = strings.Join(strings.Split(cols[3], "_")[3:], "_")
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t%s [color=red]", name))

		connections = parseAlns(cols[5])
		fileio.WriteToFileHandle(out, fmt.Sprintf("\t\"%s\" -- %s", name, connections))
		if c > 10 {
			break
		}
	}
	fileio.WriteToFileHandle(out, "}")
	exception.PanicOnErr(out.Close())
}

func parseAlns(s string) string {
	var name string
	alns := strings.Split(s, ";")
	var sb strings.Builder
	sb.WriteString("{")
	for i := 0; i < len(alns)-1; i++ {
		sb.WriteString("\"")
		name = strings.Join(strings.Split(alns[i], "_")[6:], "_")
		sb.WriteString(name)
		sb.WriteString("\"")
		sb.WriteString(" ")
	}
	sb.WriteString("}")
	return sb.String()
}
