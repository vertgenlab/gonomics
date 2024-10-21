package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

func main() {
	var cols []string

	txt := fileio.Read("/Users/sethweaver/Downloads/hsSD/regEleFam/regEleFam.txt")
	for i := range txt {
		cols = strings.Split(txt[i], "\t")
		if cols[5] == "" {
			
			continue
		}

	}
}
