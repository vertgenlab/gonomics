package dnaThreeBit

import (
)

// Append adds "b" to the end of "fragment."
// "fragment" can be nil.
func Append(fragment *ThreeBit, b ThreeBitBase) *ThreeBit { // too confusing to have Append and append (from std lib) in this package?
	var bNum uint64 = uint64(b)
	if fragment == nil {
		return &ThreeBit{Seq: []uint64{bNum << 61}, Len: 1}
	}
	var basesInLastIdx int = fragment.Len % 21
	if basesInLastIdx == 0 { // there is no more room in the slice of uint64
		fragment.Seq = append(fragment.Seq, bNum << 61) // shift the new base all the way to the left
	} else { // there is room for at least one more base in the last uint64
		var lastIdx int = len(fragment.Seq) - 1
		fragment.Seq[lastIdx] = fragment.Seq[lastIdx] | (bNum << (61 - basesInLastIdx * 3)) // bit-wise or adds the new base
	}
	fragment.Len++
	return fragment
}

// Cat appends "b" to "a."  "a" is changed to be both of them, and "b" is unchanged.
// It is quickest to have "a" be the longer sequence.
func Cat(a *ThreeBit, b *ThreeBit) { // I was tempted to call this ligate
	if b == nil {
		return
	}
	for i := 0; i < b.Len; i++ {
		a = Append(a, GetThreeBitBase(b, i))
	}
}

func Copy(a *ThreeBit) *ThreeBit {
	if a == nil {
		return nil
	}
	answer := &ThreeBit{Seq: make([]uint64, a.Len), Len: a.Len}
	copy(answer.Seq, a.Seq)
	return answer
}

