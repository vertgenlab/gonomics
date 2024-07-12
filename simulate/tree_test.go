package simulate

import (
	"math/rand"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/expandedTree"
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
		Expected:   "(((Child_2:0.014005,Child_1:0.076481)Child_6:0.057662,(Child_4:0.104908,Child_3:0.021556)Child_5:0.019682)Child_8:0.027327,Child_7:0.102769)root:0.000000;",
	},
	{NumNodes: 27,
		GammaAlpha: 1,
		GammaBeta:  50,
		SetSeed:    29,
		Expected:   "((Child_22:0.020201,((Child_8:0.002304,Child_7:0.038078)Child_10:0.000399,Child_9:0.020751)Child_21:0.024036)Child_26:0.014117,((((Child_12:0.019703,Child_11:0.000150)Child_16:0.007303,((Child_6:0.029961,(Child_4:0.006149,Child_3:0.022247)Child_5:0.003846)Child_14:0.057187,Child_13:0.078920)Child_15:0.013370)Child_18:0.000546,(Child_2:0.050898,Child_1:0.022168)Child_17:0.019087)Child_24:0.003478,(Child_20:0.016904,Child_19:0.029799)Child_23:0.036332)Child_25:0.001273)root:0.000000;",
	},
	{NumNodes: 109,
		GammaAlpha: 1.5,
		GammaBeta:  40,
		SetSeed:    31,
		Expected:   "(((Child_100:0.045150,(Child_90:0.006654,((Child_20:0.095028,Child_19:0.030417)Child_44:0.007692,Child_43:0.019950)Child_89:0.052328)Child_99:0.016959)Child_104:0.020632,(((Child_48:0.015577,Child_47:0.003022)Child_92:0.009921,(Child_84:0.050427,(Child_70:0.032997,Child_69:0.018735)Child_83:0.013334)Child_91:0.031396)Child_94:0.051754,(Child_56:0.007468,(Child_12:0.041557,Child_11:0.034014)Child_55:0.015529)Child_93:0.063935)Child_103:0.011651)Child_108:0.025138,((((Child_30:0.017368,Child_29:0.044994)Child_76:0.020903,(((Child_46:0.011248,(Child_36:0.035140,(Child_2:0.026579,Child_1:0.053884)Child_35:0.007590)Child_45:0.036926)Child_66:0.017160,Child_65:0.033690)Child_74:0.090408,Child_73:0.003538)Child_75:0.101039)Child_102:0.033161,(((Child_22:0.034958,Child_21:0.022938)Child_78:0.028343,((((Child_14:0.013471,Child_13:0.128199)Child_18:0.002494,Child_17:0.076447)Child_34:0.020413,Child_33:0.044467)Child_72:0.010106,(Child_38:0.059240,(Child_32:0.014845,(Child_6:0.017293,(Child_4:0.016716,Child_3:0.042388)Child_5:0.040786)Child_31:0.001768)Child_37:0.030235)Child_71:0.108227)Child_77:0.009140)Child_98:0.133080,(((Child_26:0.036403,Child_25:0.075220)Child_64:0.063180,Child_63:0.024564)Child_96:0.027381,((Child_50:0.015756,(Child_10:0.029784,Child_9:0.029422)Child_49:0.038286)Child_80:0.103758,(Child_68:0.004215,((Child_42:0.079636,(Child_40:0.003579,Child_39:0.054761)Child_41:0.023552)Child_62:0.013845,(Child_58:0.034503,Child_57:0.007185)Child_61:0.021921)Child_67:0.029096)Child_79:0.017513)Child_95:0.011634)Child_97:0.008438)Child_101:0.001565)Child_106:0.002624,(((Child_8:0.050998,Child_7:0.043461)Child_86:0.004088,((((Child_16:0.001782,Child_15:0.009845)Child_24:0.030671,Child_23:0.000859)Child_52:0.062528,Child_51:0.025193)Child_82:0.033179,Child_81:0.102028)Child_85:0.031076)Child_88:0.025543,(Child_60:0.040104,((Child_28:0.045015,Child_27:0.007059)Child_54:0.019888,Child_53:0.019221)Child_59:0.017869)Child_87:0.018760)Child_105:0.018000)Child_107:0.103645)root:0.000000;",
	},
}

func TestETree(t *testing.T) {
	var observed string
	for _, v := range ETreeTests {
		rand.New(rand.NewSource(v.SetSeed))
		observed = expandedTree.ToNewickString(ETree(v.NumNodes, v.GammaAlpha, v.GammaBeta))
		if strings.Compare(observed, v.Expected) != 0 {
			t.Log(observed)
			t.Errorf("Error: simulate.ETree did not produce the expected output.")
		}
	}
}
