package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

func main() {
	var j int = 0
	var k int = 0
	var bit int32
	var umis []string

	ch, _ := sam.GoReadToChan("/Users/sethweaver/Downloads/2dayBrainIueOut/lowe/possorted_genome_bam.chrSS.bam")

	for i := range ch {
		if j >= 10 {
			break
		}
		num, _, _ := sam.QueryTag(i, "xf")
		bit = num.(int32)

		if bit&8 == 8 {
			construct, _, _ := sam.QueryTag(i, "GX")
			name := construct.(string)
			umis = append(umis, name)
			k++
		}
		//j++
	}
	fmt.Println(k)
	fileio.Write("/Users/sethweaver/Downloads/2dayBrainIueOut/lowe/cellrangerUsedUMIs.txt", umis)
}
