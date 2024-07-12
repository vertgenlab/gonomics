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
		Expected:   "((Child_6:0.066447,(Child_4:0.066447,(Child_2:0.066447,Child_1:0.014881)Child_3:0.014881)Child_5:0.014881)Child_8:0.066447,Child_7:0.014881)root:0.000000;",
	},
	{NumNodes: 27,
		GammaAlpha: 1,
		GammaBeta:  50,
		SetSeed:    29,
		Expected:   "((Child_12:0.009243,(Child_10:0.009243,Child_9:0.000394)Child_11:0.000394)Child_26:0.009243,((Child_2:0.009243,Child_1:0.000394)Child_24:0.009243,(((Child_18:0.009243,((Child_6:0.009243,(Child_4:0.009243,Child_3:0.000394)Child_5:0.000394)Child_16:0.009243,(Child_8:0.009243,Child_7:0.000394)Child_15:0.000394)Child_17:0.000394)Child_20:0.009243,(Child_14:0.009243,Child_13:0.000394)Child_19:0.000394)Child_22:0.009243,Child_21:0.000394)Child_23:0.000394)Child_25:0.000394)root:0.000000;",
	},
	{NumNodes: 109,
		GammaAlpha: 1.5,
		GammaBeta:  40,
		SetSeed:    31,
		Expected:   "((((Child_98:0.022205,(Child_96:0.022205,(((Child_32:0.022205,Child_31:0.003182)Child_40:0.022205,(Child_20:0.022205,Child_19:0.003182)Child_39:0.003182)Child_44:0.022205,Child_43:0.003182)Child_95:0.003182)Child_97:0.003182)Child_100:0.022205,(((Child_82:0.022205,((Child_34:0.022205,(Child_28:0.022205,(Child_10:0.022205,Child_9:0.003182)Child_27:0.003182)Child_33:0.003182)Child_78:0.022205,((Child_36:0.022205,(Child_30:0.022205,Child_29:0.003182)Child_35:0.003182)Child_74:0.022205,Child_73:0.003182)Child_77:0.003182)Child_81:0.003182)Child_84:0.022205,(Child_72:0.022205,(Child_68:0.022205,Child_67:0.003182)Child_71:0.003182)Child_83:0.003182)Child_86:0.022205,((Child_62:0.022205,(Child_50:0.022205,Child_49:0.003182)Child_61:0.003182)Child_80:0.022205,((Child_60:0.022205,Child_59:0.003182)Child_76:0.022205,Child_75:0.003182)Child_79:0.003182)Child_85:0.003182)Child_99:0.003182)Child_102:0.022205,Child_101:0.003182)Child_108:0.022205,(((Child_92:0.022205,(((Child_8:0.022205,Child_7:0.003182)Child_58:0.022205,Child_57:0.003182)Child_90:0.022205,(((Child_12:0.022205,Child_11:0.003182)Child_48:0.022205,((Child_22:0.022205,Child_21:0.003182)Child_26:0.022205,Child_25:0.003182)Child_47:0.003182)Child_56:0.022205,(Child_46:0.022205,Child_45:0.003182)Child_55:0.003182)Child_89:0.003182)Child_91:0.003182)Child_94:0.022205,Child_93:0.003182)Child_106:0.022205,(((((Child_4:0.022205,Child_3:0.003182)Child_16:0.022205,Child_15:0.003182)Child_66:0.022205,(Child_2:0.022205,Child_1:0.003182)Child_65:0.003182)Child_88:0.022205,((((Child_18:0.022205,Child_17:0.003182)Child_52:0.022205,((Child_38:0.022205,Child_37:0.003182)Child_42:0.022205,((Child_6:0.022205,Child_5:0.003182)Child_24:0.022205,Child_23:0.003182)Child_41:0.003182)Child_51:0.003182)Child_54:0.022205,Child_53:0.003182)Child_64:0.022205,(Child_14:0.022205,Child_13:0.003182)Child_63:0.003182)Child_87:0.003182)Child_104:0.022205,(Child_70:0.022205,Child_69:0.003182)Child_103:0.003182)Child_105:0.003182)Child_107:0.003182)root:0.000000;",
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
