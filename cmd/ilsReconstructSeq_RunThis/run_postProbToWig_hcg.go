package main

import (
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/maf"
	"fmt"
)

func main() {
    fmt.Println("Hello")
	postProbFile := "/net/bmc-lab4/data/kellis/users/sarahaz/data/pdna/hcg_trails_reconstruct/post_prob.csv"
	inMafFile := "/net/bmc-lab4/data/kellis/users/rimangan/primateT2T/alignment/work/yl726/PrimateT2T_15way/outputb/chr1.maf"
	inMaf := maf.Read(inMafFile)
	out := reconstruct.PostProbToWig(postProbFile, inMaf)
	wig.Write("/net/bmc-lab4/data/kellis/users/sarahaz/data/pdna/hcg_trails_reconstruct/hcg_postProb.wig", out)
	
	fmt.Println("")
	fmt.Println("DONE")
	fmt.Println("DONE")
	fmt.Println("DONE")
}