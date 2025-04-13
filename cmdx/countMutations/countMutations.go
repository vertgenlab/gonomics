package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
)

func main() {
	var cig []cigar.Cigar
	var j, pos int

	flag.Parse()

	out := fileio.EasyCreate(flag.Arg(2))
	fileio.WriteToFileHandle(out, "readName\tmutationType\tpostion\tlength")

	in, _ := sam.GoReadToChan(flag.Arg(0))
	thresh := parse.StringToInt(flag.Arg(1))

	for i := range in {
		if i.Flag&256 == 256 || i.Flag&4 == 4 {
			continue
		}
		pos = i.GetChromStart()
		cig = i.Cigar
		for j = range cig {
			if (cig[j].Op == cigar.Deletion || cig[j].Op == cigar.Insertion) && cig[j].RunLength >= thresh {
				fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%c\t%d\t%d", i.QName, cig[j].Op, pos, cig[j].RunLength))
			}
			if cigar.ConsumesReference(cig[j].Op) {
				pos += cig[j].RunLength
			}
		}
	}
	exception.PanicOnErr(out.Close())
}
