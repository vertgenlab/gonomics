package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/tree_newick"
	"math/rand"
	"strings"
	"time"
)

// makes random gene with start and stop codon, must be length multiple of 3
func Rand_gene(name string, length int) {
	if length%3 != 0 {
		fmt.Print("length must be divisible by three")
	} else {
		seq := "ATG"
		rand_length := length - 6
		for i := 0; i < rand_length; i++ {
			rt := rand.NewSource(time.Now().UnixNano())
			rn := rand.New(rt)
			r := rn.Intn(100)

			switch {
			case r < 25:
				seq = seq + "A"
			case r < 50:
				seq = seq + "C"
			case r < 75:
				seq = seq + "T"
			default:
				seq = seq + "G"

			}
		}

		rts := rand.NewSource(time.Now().UnixNano())
		rns := rand.New(rts)
		rs := rns.Intn(100)
		switch {
		case rs < 100/3:
			seq = seq + "TAG"
		case rs < 100/3*2:
			seq = seq + "TGA"
		default:
			seq = seq + "TAA"
		}

		sb := dna.StringToBases(seq)

		record := fasta.Fasta{name, sb}
		var answer []*fasta.Fasta
		answer = append(answer, &record)
		filename := name + ".fasta"
		fasta.Write(filename, answer)
	}
}

//final function to run to simulate based off of the random gene and the tree
func Simulate(rand_seq_filename string, tree_output_filename string, root *tree_newick.NTree) {
	rand1 := fasta.Read(rand_seq_filename)
	root.Fasta = rand1[0]
	seq := Get_seq(rand1)
	fasta.Write(tree_output_filename, Tree_print(root, seq))
}

func Translate(codon []string) aa {
	cod := strings.Join(codon, "")
	var am aa
	if len(cod) != 3 {
		fmt.Print("codon must have length of 3 bases")
	} else {
		am = m[cod]
	}
	return am
}

type aa int

// BLOSUM matrix for amino acid switching probabilities
var BLOSUM = [][]float64{[]float64{0.288590604, 0.03087248322, 0.03087248322, 0.02953020134, 0.02147651007, 0.0255033557, 0.04026845638, 0.07785234899, 0.01476510067, 0.04295302013, 0.05906040268, 0.04429530201, 0.01744966443, 0.02147651007, 0.02953020134, 0.08456375839, 0.04966442953, 0.005369127517, 0.01744966443, 0.06845637584, 0.0},
	[]float64{0.04457364341, 0.3449612403, 0.03875968992, 0.03100775194, 0.007751937984, 0.0484496124, 0.0523255814, 0.03294573643, 0.02325581395, 0.02325581395, 0.04651162791, 0.1201550388, 0.01550387597, 0.01744186047, 0.01937984496, 0.04457364341, 0.03488372093, 0.005813953488, 0.01744186047, 0.03100775194, 0.0},
	[]float64{0.05122494432, 0.04454342984, 0.3140311804, 0.08240534521, 0.008908685969, 0.03340757238, 0.04899777283, 0.06458797327, 0.03118040089, 0.02227171492, 0.03118040089, 0.05345211581, 0.01113585746, 0.01781737194, 0.02004454343, 0.06904231626, 0.04899777283, 0.004454342984, 0.01559020045, 0.02672605791, 0.0},
	[]float64{0.04104477612, 0.02985074627, 0.06902985075, 0.3973880597, 0.007462686567, 0.02985074627, 0.09141791045, 0.04664179104, 0.01865671642, 0.0223880597, 0.02798507463, 0.0447761194, 0.009328358209, 0.01492537313, 0.0223880597, 0.05223880597, 0.03544776119, 0.003731343284, 0.01119402985, 0.02425373134, 0.0},
	[]float64{0.06504065041, 0.0162601626, 0.0162601626, 0.0162601626, 0.4837398374, 0.01219512195, 0.0162601626, 0.0325203252, 0.008130081301, 0.04471544715, 0.06504065041, 0.02032520325, 0.0162601626, 0.02032520325, 0.0162601626, 0.0406504065, 0.03658536585, 0.00406504065, 0.01219512195, 0.05691056911, 0.0},
	[]float64{0.05588235294, 0.07352941176, 0.04411764706, 0.04705882353, 0.008823529412, 0.2147058824, 0.1029411765, 0.04117647059, 0.02941176471, 0.02647058824, 0.04705882353, 0.09117647059, 0.02058823529, 0.01470588235, 0.02352941176, 0.05588235294, 0.04117647059, 0.005882352941, 0.02058823529, 0.03529411765, 0.0},
	[]float64{0.05524861878, 0.04972375691, 0.04051565378, 0.09023941068, 0.007366482505, 0.06445672192, 0.2965009208, 0.0349907919, 0.02578268877, 0.02209944751, 0.03683241252, 0.07550644567, 0.01289134438, 0.01657458564, 0.02578268877, 0.05524861878, 0.03683241252, 0.005524861878, 0.01657458564, 0.03130755064, 0.0},
	[]float64{0.07827260459, 0.02294197031, 0.03913630229, 0.03373819163, 0.01079622132, 0.01889338731, 0.02564102564, 0.5101214575, 0.01349527665, 0.01889338731, 0.02834008097, 0.03373819163, 0.009446693657, 0.01619433198, 0.01889338731, 0.05128205128, 0.02968960864, 0.005398110661, 0.01079622132, 0.02429149798, 0.0},
	[]float64{0.04198473282, 0.04580152672, 0.0534351145, 0.03816793893, 0.007633587786, 0.03816793893, 0.0534351145, 0.03816793893, 0.3549618321, 0.02290076336, 0.03816793893, 0.04580152672, 0.01526717557, 0.03053435115, 0.01908396947, 0.04198473282, 0.02671755725, 0.007633587786, 0.0572519084, 0.02290076336, 0.0},
	[]float64{0.0471281296, 0.0176730486, 0.0147275405, 0.0176730486, 0.01620029455, 0.01325478645, 0.0176730486, 0.0206185567, 0.0088365243, 0.2709867452, 0.1678939617, 0.0235640648, 0.03681885125, 0.0441826215, 0.0147275405, 0.02503681885, 0.03976435935, 0.0058910162, 0.0206185567, 0.176730486, 0.0},
	[]float64{0.04453441296, 0.02429149798, 0.01417004049, 0.01518218623, 0.01619433198, 0.01619433198, 0.02024291498, 0.02125506073, 0.01012145749, 0.1153846154, 0.3755060729, 0.02530364372, 0.0495951417, 0.05465587045, 0.01417004049, 0.02429149798, 0.03340080972, 0.007085020243, 0.02226720648, 0.09615384615, 0.0},
	[]float64{0.05699481865, 0.1070811744, 0.0414507772, 0.0414507772, 0.008635578584, 0.05354058722, 0.07081174439, 0.04317789292, 0.0207253886, 0.02763385147, 0.04317789292, 0.2780656304, 0.01554404145, 0.01554404145, 0.02763385147, 0.05354058722, 0.03972366149, 0.00518134715, 0.01727115717, 0.03281519862, 0.0},
	[]float64{0.05220883534, 0.03212851406, 0.02008032129, 0.02008032129, 0.01606425703, 0.0281124498, 0.0281124498, 0.0281124498, 0.01606425703, 0.1004016064, 0.1967871486, 0.03614457831, 0.1606425703, 0.04819277108, 0.01606425703, 0.03614457831, 0.04016064257, 0.008032128514, 0.02409638554, 0.09236947791, 0.0},
	[]float64{0.03382663848, 0.01902748414, 0.01691331924, 0.01691331924, 0.01057082452, 0.01057082452, 0.01902748414, 0.02536997886, 0.01691331924, 0.06342494715, 0.1141649049, 0.01902748414, 0.02536997886, 0.3868921776, 0.01057082452, 0.02536997886, 0.02536997886, 0.01691331924, 0.088794926, 0.05496828753, 0.0},
	[]float64{0.05684754522, 0.02583979328, 0.02325581395, 0.03100775194, 0.01033591731, 0.02067183463, 0.03617571059, 0.03617571059, 0.01291989664, 0.02583979328, 0.03617571059, 0.04134366925, 0.01033591731, 0.01291989664, 0.4935400517, 0.04392764858, 0.03617571059, 0.002583979328, 0.01291989664, 0.03100775194, 0.0},
	[]float64{0.109947644, 0.04013961606, 0.05410122164, 0.04886561955, 0.01745200698, 0.03315881326, 0.05235602094, 0.06631762653, 0.01919720768, 0.02966841187, 0.04188481675, 0.05410122164, 0.01570680628, 0.02094240838, 0.02966841187, 0.219895288, 0.08202443281, 0.005235602094, 0.01745200698, 0.04188481675, 0.0},
	[]float64{0.07297830375, 0.03550295858, 0.04339250493, 0.03747534517, 0.01775147929, 0.02761341223, 0.03944773176, 0.04339250493, 0.01380670611, 0.05325443787, 0.0650887574, 0.04536489152, 0.01972386588, 0.02366863905, 0.02761341223, 0.09270216963, 0.2465483235, 0.005917159763, 0.01775147929, 0.07100591716, 0.0},
	[]float64{0.0303030303, 0.02272727273, 0.01515151515, 0.01515151515, 0.007575757576, 0.01515151515, 0.02272727273, 0.0303030303, 0.01515151515, 0.0303030303, 0.05303030303, 0.02272727273, 0.01515151515, 0.06060606061, 0.007575757576, 0.02272727273, 0.02272727273, 0.4924242424, 0.06818181818, 0.0303030303, 0.0},
	[]float64{0.04049844237, 0.02803738318, 0.02180685358, 0.01869158879, 0.009345794393, 0.02180685358, 0.02803738318, 0.02492211838, 0.04672897196, 0.04361370717, 0.06853582555, 0.03115264798, 0.01869158879, 0.1308411215, 0.01557632399, 0.03115264798, 0.02803738318, 0.02803738318, 0.3177570093, 0.04672897196, 0.0},
	[]float64{0.06995884774, 0.0219478738, 0.01646090535, 0.01783264746, 0.01920438957, 0.01646090535, 0.02331961591, 0.02469135802, 0.008230452675, 0.1646090535, 0.1303155007, 0.02606310014, 0.03155006859, 0.03566529492, 0.01646090535, 0.0329218107, 0.04938271605, 0.00548696845, 0.02057613169, 0.268861454, 0.0},
	[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}

//amino acids
const (
	ala  aa = 0
	arg  aa = 1
	asn  aa = 2
	asp  aa = 3
	cys  aa = 4
	gln  aa = 5
	glu  aa = 6
	gly  aa = 7
	his  aa = 8
	ile  aa = 9
	leu  aa = 10
	lys  aa = 11
	met  aa = 12
	phe  aa = 13
	pro  aa = 14
	ser  aa = 15
	thr  aa = 16
	trp  aa = 17
	tyr  aa = 18
	val  aa = 19
	stop aa = 20
)

//translation map
var m = map[string]aa{
	"TGA": aa(20), "TAA": aa(20), "TAG": aa(20),
	"GTA": aa(19), "GTC": aa(19), "GTG": aa(19), "GTT": aa(19),
	"TAT": aa(18), "TAC": aa(18),
	"TGG": aa(17),
	"ACA": aa(16), "ACG": aa(16), "ACT": aa(16), "ACC": aa(16),
	"TCA": aa(15), "TCC": aa(15), "TCG": aa(15), "TCT": aa(15), "AGT": aa(15), "AGC": aa(15),
	"CCC": aa(14), "CCT": aa(14), "CCA": aa(14), "CCG": aa(14),
	"TTT": aa(13), "TTC": aa(13),
	"ATG": aa(12),
	"AAA": aa(11), "AAG": aa(11),
	"TTA": aa(10), "TTG": aa(10), "CTC": aa(10), "CTG": aa(10), "CTA": aa(10), "CTT": aa(10),
	"ATT": aa(9), "ATC": aa(9), "ATA": aa(9),
	"CAT": aa(8), "CAC": aa(8),
	"GGG": aa(7), "GGA": aa(7), "GGT": aa(7), "GGC": aa(7),
	"GAA": aa(6), "GAG": aa(6),
	"CAA": aa(5), "CAG": aa(5),
	"TGT": aa(4), "TGC": aa(4),
	"GAT": aa(3), "GAC": aa(3),
	"AAT": aa(2), "AAC": aa(2),
	"AGA": aa(1), "AGG": aa(1), "CGC": aa(1), "CGG": aa(1), "CGA": aa(1), "CGT": aa(1),
	"GCA": aa(0), "GCG": aa(0), "GCT": aa(0), "GCC": aa(0)}

//get sequence of fastas
func Get_seq(record []*fasta.Fasta) []string {
	var sequence []string

	for _, rec := range record {
		for i := 0; i < len(rec.Seq); i += len(rec.Seq) {
			s := dna.BasesToString(rec.Seq[i:])
			sequence = strings.Split(s, "")
		}
	}
	return sequence
}

//choose base based off of random seed
func Base_choice() string {

	rt := rand.NewSource(time.Now().UnixNano())
	rn := rand.New(rt)
	r := rn.Intn(100)
	var base string
	switch {
	case r < 25:
		base = "A"
	case r < 50:
		base = "C"
	case r < 75:
		base = "T"
	default:
		base = "G"
	}
	return base
}

//mutate base given random seeded float
func Mutate_base(b string, distance_percent float64) string {

	rt := rand.NewSource(time.Now().UnixNano())
	rn := rand.New(rt)
	r := rn.Float64()

	var base string
	switch {
	case distance_percent == 0:
		base = b
	case r < distance_percent:
		base = Base_choice()
	default:
		base = b
	}

	return base

}

//mutate sequence taking BLOSUM probabilites and gene structure into account
func Mutate_seq(seq []string, distance_percent float64) []string {
	c := make([]string, len(seq))
	copy(c, seq)
	if len(seq)%3 != 0 {
		fmt.Print("sequence length must be divisible by three")
	} else {
		rt := rand.NewSource(time.Now().UnixNano())
		rn := rand.New(rt)
		r := rn.Float64()

		var base string
		var base_new string
		var codon []string
		var codon_new []string
		var am_ac aa
		var am_ac_new aa
		codons := len(seq) / 3
		for i := 0; i < codons; i++ {
			codon = []string{seq[i*3], seq[i*3+1], seq[i*3+2]}
			for j := 0; j < 3; j++ {
				base = c[i*3+j]
				switch {
				case i == 0:
					base_new = base
				case i == codons:
					base_new = base
				default:
					base_new = Mutate_base(base, distance_percent)
				}
				codon_new = []string{seq[i*3], seq[i*3+1], seq[i*3+2]}
				codon_new[j] = base_new
				am_ac = Translate(codon)

				am_ac_new = Translate(codon_new)
				if am_ac_new != am_ac {
					prob := BLOSUM[am_ac][am_ac_new]
					if r < prob {

						c[i*3+j] = base_new
					} else {
						c[i*3+j] = base

					}
				} else {
					c[i*3+j] = base_new
				}
			}

		}

	}
	return c
}

func SeqtoFformat(sequence []string, new_name string) fasta.Fasta {
	seq := strings.Join(sequence, "")
	var answer fasta.Fasta
	for i := 0; i < len(seq); i++ {
		s := dna.StringToBases(seq)
		answer = fasta.Fasta{new_name, s}
	}
	return answer
}

//make fastas based off of node and random DNA_sequence
func Tree_print(node *tree_newick.NTree, DNA_seq []string) []*fasta.Fasta {
	var fasta_final []*fasta.Fasta
	var seq []string

	length := float64(node.BranchLength)

	seq = Mutate_seq(DNA_seq, length)

	s := SeqtoFformat(seq, node.Name)
	fasta_final = append(fasta_final, &s)
	if node.Left != nil && node.Right != nil {
		b := Tree_print(node.Right, seq)
		fasta_final = append(fasta_final, b...)
		a := Tree_print(node.Left, seq)
		fasta_final = append(fasta_final, a...)
	}
	return fasta_final
}

//remove ancestors from the fasta file for reconstruction
func Remove_No_Names(filename string) {
	fastas := fasta.Read(filename)
	var fastas_new []*fasta.Fasta
	for i := 0; i < len(fastas); i++ {
		if fastas[i].Name != "" {
			fastas_new = append(fastas_new, fastas[i])
		}

	}
	filename_new := "descendents_" + filename
	fasta.Write(filename_new, fastas_new)
}

func Get_leaf(node *tree_newick.NTree) []*tree_newick.NTree {
	var leaf []*tree_newick.NTree
	if node.Left != nil && node.Right != nil {
		a := Get_leaf(node.Left)
		b := Get_leaf(node.Right)
		leaf = append(leaf, a...)
		leaf = append(leaf, b...)
	}
	if node.Left == nil && node.Right == nil {
		leaf = append(leaf, node)
	}
	return leaf
}

func Remove_ancestors(filename string, tree *tree_newick.NTree) {
	fastas := fasta.Read(filename)
	var fastas_new []*fasta.Fasta

	leaf := Get_leaf(tree)
	for i := 0; i < len(fastas); i++ {
		for j := 0; j < len(leaf); j++ {
			if fastas[i].Name == leaf[j].Name {
				fastas_new = append(fastas_new, fastas[i])
			}
		}
	}
	filename_new := "descendents_" + filename
	fasta.Write(filename_new, fastas_new)
}
