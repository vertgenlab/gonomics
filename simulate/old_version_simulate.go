package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sophie/reconstruct"
	"math/rand"
	"strings"
	"time"
)

func Rand_seq(name string, length int) {
	seq := ""
	for i := 0; i < length; i++ {
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
	sb, err := dna.StringToBases(seq)
	fmt.Print(err)

	record := fasta.Fasta{name, sb}
	var answer []fasta.Fasta
	answer = append(answer, record)
	filename := name + ".fasta"
	fasta.Write(filename, answer)

}
func Get_seq(record []fasta.Fasta) []string {
	var sequence []string

	for _, rec := range record {
		for i := 0; i < len(rec.Seq); i += len(rec.Seq) {
			s := dna.BasesToString(rec.Seq[i:])
			sequence = strings.Split(s, "")
		}
	}

	return sequence
}

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

func Mutate_seq(seq []string, distance_percent float64) []string {

	c := make([]string, len(seq))
	copy(c, seq)

	var base string
	var base_new string

	for i := 0; i < len(seq); i++ {
		base = c[i]
		base_new = Mutate_base(base, distance_percent)
		c[i] = base_new

	}
	// fmt.Print(c, "\n")
	return c

}

func SeqtoFformat(sequence []string, new_name string) fasta.Fasta {
	seq := strings.Join(sequence, "")
	var new fasta.Fasta
	for i := 0; i < len(seq); i++ {
		s, err := dna.StringToBases(seq)

		if err == nil {
			new = fasta.Fasta{new_name, s}

		}
	}
	return new
}

func Tree_print(node *reconstruct.NTree, DNA_seq []string) []fasta.Fasta {
	var fasta_final []fasta.Fasta
	var seq []string

	length := float64(node.BranchLength)

	seq = Mutate_seq(DNA_seq, length)
	// fmt.Print(seq, "\n")

	s := SeqtoFformat(seq, node.Name)
	fasta_final = append(fasta_final, s)

	if node.Left != nil && node.Right != nil {
		b := Tree_print(node.Right, seq)
		fasta_final = append(fasta_final, b...)
		a := Tree_print(node.Left, seq)
		fasta_final = append(fasta_final, a...)
	}

	return fasta_final

}

func Remove_ancestors(filename string) {
	fastas, err := fasta.Read(filename)
	if err != nil {
	}
	var fastas_new []fasta.Fasta
	for i := 0; i < len(fastas); i++ {
		if fastas[i].Name != "" {
			fastas_new = append(fastas_new, fastas[i])
		}

	}
	filename_new := "descendents_" + filename
	fasta.Write(filename_new, fastas_new)
}

func Simulate(rand_seq_filename string, tree_output_filename string, root *reconstruct.NTree) {
	rand1, err := fasta.Read(rand_seq_filename)
	fmt.Print(err)
	root.Fasta = &rand1[0]
	seq := Get_seq(rand1)
	fasta.Write(tree_output_filename, Tree_print(root, seq))
}
