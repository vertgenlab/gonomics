package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

func main() {
	var k int = 0
	var bit int32
	var umis []string

	//ch, _ := sam.GoReadToChan("/Users/sethweaver/Downloads/2dayBrainIueOut/lowe/possorted_genome_bam.chrSS.bam") //lowe library analysis
	//ch, _ := sam.GoReadToChan("/Users/sethweaver/Downloads/2dayBrainIueOut/silver/silver_possorted_genome_bam.chrSS.bam") //silver library analysis
	ch, _ := sam.GoReadToChan("/Users/sethweaver/Downloads/2dayBrainIueOut/haqer/haqer_possorted_genome_bam.chrSS.bam") //haqer reanalysis

	for i := range ch {
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags
		bit = num.(int32)

		if bit&8 == 8 { // bit 8 is the flag for UMI used in final count.
			construct, _, _ := sam.QueryTag(i, "GX") // gene associated with read
			name := construct.(string)
			umis = append(umis, name)
			k++
		}
	}
	fmt.Println(k)
	fileio.Write("/Users/sethweaver/Downloads/2dayBrainIueOut/haqer/cellrangerUsedUMIs.txt", umis)
}
