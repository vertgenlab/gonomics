package giraf

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

type Giraf struct {
	QName     string
	QStart    int
	QEnd      int
	PosStrand bool
	Path      *Path
	Aln       *cigar.Cigar // current cigar will need to be expanded
	AlnScore  int
	MapQ      uint8
	Seq       []dna.Base // dnaTwoBit?
	Qual      []uint8
	Notes     []Note // Similar to sam, this is should be a list of notes.
	// Each note should be of the form TAG:TYPE:VALUE
	// TAG is two characters
	// TYPE is a single character
	// VALUE will be stored as a string and can then be de-coded based on type
	// An example would be "BZ:4000
}

type Path struct {
	TStart int      // The path starts on the TStart base (0-based, closed) of Nodes[0]
	Nodes  []uint32 // The node Id/Index of all the nodes in the path
	TEnd   int      // The path ends on the TEnd base (0-based, open) of Nodes[len(Nodes)-1]
}

type Note struct {
	Tag   string
	Type  rune
	Value string
}
