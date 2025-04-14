package main

import (
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/wig"
)

func main() {
    fmt.Println("Hello")
	postProbFile := "testdata/hcg_post_prob.csv"
	inMaf := "/net/bmc-lab4/data/kellis/users/rimangan/primateT2T/alignment/work/yl726/PrimateT2T_15way/outputb/chr1.maf"
	out = PostProbToWig(postProbFile, inMaf)
	wig.Write("/net/bmc-lab4/data/kellis/users/sarahaz/data/pdna/hcg_trails_reconstruct/hcg_postProb.wig", out)
	
	fmt.Println("")
	fmt.Println("DONE")
	fmt.Println("DONE")
	fmt.Println("DONE")
}