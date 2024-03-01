package genomeGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

const (
	left  byte = 0
	right byte = 1
)

// getTargetBases retrieves a slice of bases from the specified direction, either 'left' or 'right'.
func getTargetBases(n *Node, extension, position int, seq []dna.Base, ans []dna.Base, direction byte) []dna.Base {
	var basesToTake int
	if direction == left {
		// Calculate the number of bases we can take from the left.
		basesToTake = numbers.Min(len(seq)+position, extension) - len(seq)
	} else if direction == right {
		// Calculate the number of bases we can take from the right.
		basesToTake = numbers.Min(len(seq)+len(n.Seq)-position, extension) - len(seq)
	}

	// Ensure 'ans' has enough capacity to avoid reallocation.
	requiredLen := len(ans) + len(seq) + basesToTake
	if cap(ans) < requiredLen {
		// Allocate a new slice with enough capacity.
		newAns := make([]dna.Base, requiredLen, requiredLen+extension)
		copy(newAns, ans)
		ans = newAns
	} else {
		// Use available capacity in 'ans'.
		ans = ans[:requiredLen]
	}

	if direction == left {
		// Copy the bases from the left into 'ans'.
		copy(ans, n.Seq[position-basesToTake:position])
		copy(ans[basesToTake:], seq)
	} else if direction == right {
		// Copy the bases from the right into 'ans'.
		copy(ans, seq)
		copy(ans[len(seq):], n.Seq[position:position+basesToTake])
	}

	return ans
}
