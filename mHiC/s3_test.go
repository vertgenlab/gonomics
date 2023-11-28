package mHiC

import (
	"testing"
)

func TestGetPrefix(t *testing.T) {
	in := "foo.bar.bar"
	exp := "foo.bar"
	out := getPrefix(in)
	if exp != out {
		t.Errorf("Exp: %s\tOut: %s\n", exp, out)
	}
}
