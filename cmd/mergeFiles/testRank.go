package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

func main() {
	var c int
	var a, b sam.Sam
	//samFile, _ := sam.Read("testdata/small.sam")
	//samFile, _ := sam.Read("testdata/small.sortName.bam")
	samFile, _ := sam.GoReadToChan("testdata/test1.sort.sam")

	for i := 0; i < 1000000000; i++ {
		a = <-samFile
		b = <-samFile
		ans := strings.Compare(a.QName, b.QName)
		if ans == 1 {
			//fmt.Println(a.QName, b.QName)
			c++
		}
	}
	fmt.Println(c)
}
