package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/vcf"
)

type test struct {
	i int
	s string
	b bool
	p *test
}

func main() {
	/*chroms := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
		"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
		"chr22", "chrX", "chrY"}

	for i := range chroms {
		fmt.Printf("lastz ../hs1.byChrom/%s.hs1.fa --self --output=withMask/axtTest/selfAlign/%s_%s.hs1.lav --allocate:traceback=1G --scores=human_chimp_v2.mat O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0\n", chroms[i], chroms[i], chroms[i])
		//for j := i + 1; j < len(chroms); j++ {
		//	fmt.Printf("lastz ../hs1.byChrom/%s.hs1.fa ../hs1.byChrom/%s.hs1.fa --format=axt --output=withMask/axtTest/%s_%s.hs1.axt --allocate:traceback=1G --scores=human_chimp_v2.mat O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0\n", chroms[i], chroms[j], chroms[i], chroms[j])
		//}
	}


		var cols []string
		txt := fileio.Read("/Users/sethweaver/Downloads/hsSD/netSdDiscovery/pt.chromNames.txt")
		for i := range txt {
			cols = strings.Split(txt[i], "\t")
			fmt.Printf("sed -i 's/%s/%s/g hs1.chainSynGCA_028858775.2.bed \n", cols[0], cols[1])
		}


	*/

	var snps []vcf.Vcf
	var size int
	var i int

	flag.Parse()

	fmt.Println("length\tpercID")
	xt, _ := axt.GoReadToChan(flag.Arg(0))

	for rec := range xt {
		i++
		snps = axt.ToVcf(rec)
		size = rec.GetChromEnd() - rec.GetChromStart()
		fmt.Printf("%d\t%f\n", size, 100-(float64(len(snps))/float64(size))*100)
		if i > 10 {
			break
		}
	}
}
