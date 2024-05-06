package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/expandedTree"
	"math/rand"
	"strings"
	"testing"
)

var ETreeTests = []struct {
	NumNodes   int
	GammaAlpha float64
	GammaBeta  float64
	SetSeed    int64
	Expected   string
}{
	{NumNodes: 9,
		GammaAlpha: 2,
		GammaBeta:  20,
		SetSeed:    14,
		Expected:   "((Child_6:0.034954,(Child_4:0.139110,Child_3:0.086596)Child_5:0.341105)Child_8:0.045311,(Child_2:0.166524,Child_1:0.050019)Child_7:0.123820)root:0.000000;",
	},
	{NumNodes: 27,
		GammaAlpha: 1,
		GammaBeta:  50,
		SetSeed:    29,
		Expected:   "(((Child_22:0.003211,(Child_20:0.002872,Child_19:0.042758)Child_21:0.006644)Child_24:0.037514,(((Child_10:0.013719,((Child_6:0.061806,Child_5:0.013832)Child_8:0.018914,Child_7:0.007295)Child_9:0.002391)Child_14:0.014917,((Child_2:0.051289,Child_1:0.006884)Child_4:0.015913,Child_3:0.021120)Child_13:0.014485)Child_18:0.006290,Child_17:0.008466)Child_23:0.023560)Child_26:0.036208,(Child_16:0.040451,(Child_12:0.039931,Child_11:0.003423)Child_15:0.005260)Child_25:0.003926)root:0.000000;",
	},
	{NumNodes: 109,
		GammaAlpha: 1.5,
		GammaBeta:  40,
		SetSeed:    31,
		Expected:   "((((Child_44:0.006517,Child_43:0.109418)Child_102:0.022691,Child_101:0.069907)Child_104:0.017401,(((Child_88:0.058848,Child_87:0.068688)Child_98:0.039568,((Child_48:0.004348,Child_47:0.046161)Child_66:0.085396,(Child_64:0.034291,(Child_52:0.015439,(Child_16:0.021705,(Child_4:0.004945,Child_3:0.057334)Child_15:0.053088)Child_51:0.116166)Child_63:0.005553)Child_65:0.027644)Child_97:0.025311)Child_100:0.089792,((Child_40:0.025914,Child_39:0.039614)Child_46:0.043156,(Child_34:0.044620,Child_33:0.029645)Child_45:0.023004)Child_99:0.020287)Child_103:0.064126)Child_108:0.089536,(((Child_86:0.033879,((Child_58:0.039604,(Child_42:0.045386,Child_41:0.093776)Child_57:0.027673)Child_72:0.026161,((Child_6:0.032276,Child_5:0.031596)Child_30:0.004923,Child_29:0.012088)Child_71:0.028735)Child_85:0.008908)Child_90:0.021336,((Child_62:0.036533,Child_61:0.033974)Child_84:0.054544,Child_83:0.014184)Child_89:0.025722)Child_106:0.015073,((Child_94:0.005651,(((((Child_8:0.020075,Child_7:0.021594)Child_24:0.020644,Child_23:0.023004)Child_74:0.053036,(Child_68:0.028747,Child_67:0.010411)Child_73:0.011241)Child_80:0.084856,Child_79:0.027478)Child_92:0.133742,(((Child_54:0.054877,(Child_50:0.035462,Child_49:0.002643)Child_53:0.019713)Child_78:0.047785,(((Child_14:0.035424,Child_13:0.005973)Child_38:0.058865,(((Child_28:0.015828,Child_27:0.015684)Child_32:0.003373,Child_31:0.084926)Child_36:0.060556,(Child_26:0.029531,((Child_2:0.059249,Child_1:0.039055)Child_18:0.026767,Child_17:0.004992)Child_25:0.124374)Child_35:0.031473)Child_37:0.028492)Child_76:0.012725,(Child_20:0.018221,Child_19:0.039200)Child_75:0.013974)Child_77:0.071601)Child_82:0.077139,(Child_22:0.016856,Child_21:0.004978)Child_81:0.129835)Child_91:0.041397)Child_93:0.080565)Child_96:0.041267,((Child_56:0.010595,Child_55:0.011583)Child_70:0.122753,((Child_12:0.018362,(Child_10:0.088531,Child_9:0.041212)Child_11:0.108323)Child_60:0.007330,Child_59:0.038211)Child_69:0.007207)Child_95:0.031188)Child_105:0.005361)Child_107:0.035403)root:0.000000;",
	},
}

func TestETree(t *testing.T) {
	var observed string
	for _, v := range ETreeTests {
		rand.Seed(v.SetSeed)
		observed = expandedTree.ToNewickString(ETree(v.NumNodes, v.GammaAlpha, v.GammaBeta))
		if strings.Compare(observed, v.Expected) != 0 {
			fmt.Println(observed)
			t.Errorf("Error: simulate.ETree did not produce the expected output.")
		}
	}
}
