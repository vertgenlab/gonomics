package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
)

func main() {
	var target, ans, tmpAns []interval.Interval
	synNet := bed.Read("/Users/sethweaver/Downloads/hsSD/netSdDiscovery/hs1.chainSynGCA_028858775.2.bed")
	for i := range synNet {
		target = append(target, synNet[i])
	}
	tree := interval.BuildTree(target)

	for i := range synNet {
		tmpAns = interval.Query(tree, synNet[i], "di")
		if len(tmpAns) > 0 {
			ans = append(ans, synNet[i])
		}
	}

	out := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/netSdDiscovery/2ndLvlAln.humanChimp.bed")
	for i := range ans {
		bed.WriteToFileHandle(out, ans[i].(bed.Bed))
	}
	exception.PanicOnErr(out.Close())
}
