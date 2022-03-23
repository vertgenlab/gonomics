package common

import (
	"golang.org/x/exp/constraints"
	"testing"
)

func TestGenericWorkflowSupport(t *testing.T) {
	generic(1, 2)
}

func generic[E constraints.Integer](a, b E) bool {
	return a < b
}
