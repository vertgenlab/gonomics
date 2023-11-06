package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

func mergeFiles(sam1, sam2, outFile string) {
	var currF1 sam.Sam
	var ans int

	out := fileio.EasyCreate(outFile)

	f1, _ := sam.GoReadToChan(sam1)
	f2, head := sam.GoReadToChan(sam2)

	sam.WriteHeaderToFileHandle(out, head)

	currF1 = <-f1
	for i := range f2 {
		ans = strings.Compare(i.QName, currF1.QName)
		switch ans {
		case 0:
			sam.WriteToFileHandle(out, currF1)
		case 1:
			for ans == 1 {
				currF1 = <-f1
				ans = strings.Compare(i.QName, currF1.QName)
				switch ans {
				case 0:
					sam.WriteToFileHandle(out, currF1)
				case -1:
					sam.WriteToFileHandle(out, i)
				}
			}
		case -1:
			sam.WriteToFileHandle(out, i)
		}
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func main() {
	flag.Parse()
	mergeFiles(flag.Arg(0), flag.Arg(1), flag.Arg(2))
}
