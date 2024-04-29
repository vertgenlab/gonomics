package dnaTwoBit

// Append adds a TwoBitBase to the end of a TwoBit sequence.
func Append(fragment *TwoBit, b TwoBitBase) *TwoBit {
	if fragment == nil {
		return &TwoBit{Seq: []uint64{uint64(b) << 62}, Len: 1}
	}
	var bNum uint64 = uint64(b)
	basesInLastIdx := fragment.Len % 32
	if basesInLastIdx == 0 {
		fragment.Seq = append(fragment.Seq, bNum<<62)
	} else {
		lastIdx := len(fragment.Seq) - 1
		fragment.Seq[lastIdx] |= (bNum << (62 - 2*basesInLastIdx))
	}
	fragment.Len++
	return fragment
}

// Cat concatenates two TwoBit sequences, modifying the first sequence to include both.
func Cat(a *TwoBit, b *TwoBit) { // I was tempted to call this ligate
	if b == nil {
		return
	}
	for i := 0; i < b.Len; i++ {
		a = Append(a, GetTwoBitBase(b, uint(i)))
	}
}

// Copy returns a duplicate of the TwoBit passed in.
func Copy(a *TwoBit) *TwoBit {
	if a == nil {
		return nil
	}
	copySeq := make([]uint64, len(a.Seq))
	copy(copySeq, a.Seq)
	return &TwoBit{Seq: copySeq, Len: a.Len}
}
