package tree

import (
	"testing"
)

var treeTests = []struct {
	treeText string
}{
	{"(human,chimp)ancestor;"},
	{"((human,chimp),rhesus);"},
	{"(((human,chimp),(mouse,rat)),dog);"},
	{"((platypus,(opossum,((((rat,mouse),rabbit),human),dog))),(lizard,bird));"},
	{"(((((((((((((((((((hg19:0.006653,panTro2:0.006643):0.002697,gorGor1:0.008759):0.009628,ponAbe2:0.018154):0.040041,rheMac2:0.008857):0.002523,papHam1:0.008659):0.045126,calJac1:0.066538):0.057197,tarSyr1:0.138123):0.010973,(micMur1:0.092952,otoGar1:0.129733):0.035547):0.015429,tupBel1:0.186655):0.004774,(((((mm9:0.084491,rn4:0.091658):0.198478,dipOrd1:0.212161):0.022993,cavPor3:0.225939):0.009931,speTri1:0.148899):0.025744,(oryCun2:0.114209,ochPri2:0.201057):0.101960):0.015281):0.020758,(((vicPac1:0.107324,(turTru1:0.064649,bosTau4:0.123637):0.025222):0.040511,((equCab2:0.109425,(felCat3:0.098713,canFam2:0.102642):0.050039):0.006102,(myoLuc1:0.142459,pteVam1:0.113411):0.033834):0.004353):0.011465,(eriEur1:0.222309,sorAra1:0.269795):0.056590):0.021319):0.023750,(((loxAfr3:0.082175,proCap1:0.155554):0.026780,echTel1:0.246288):0.050161,(dasNov2:0.116597,choHof1:0.096282):0.053313):0.006291):0.401417,macEug1:0.133639):0.002591,monDom5:0.150986):0.200139,ornAna1:0.459714):0.119246,((galGal3:0.163701,taeGut1:0.172622):0.200929,anoCar1:0.487683):0.101155):0.183043,xenTro2:0.837839):0.323351,(((tetNig2:0.224064,fr2:0.202940):0.192554,(gasAcu1:0.313558,oryLat2:0.480838):0.062205):0.323393,danRer6:0.730241):0.157048):0.524881,petMar1:0.524881);"},
}

func TestStringConv(t *testing.T) {
	for _, test := range treeTests {
		tree, err := parseNewick(test.treeText)
		if err != nil {
			t.Error(err)
		}
		answer := ToString(tree)
		if test.treeText != answer {
			t.Errorf("Error: after converting text to a tree and back:\n%s\nbecame\n%s\n", test.treeText, answer)
		}
	}
}
