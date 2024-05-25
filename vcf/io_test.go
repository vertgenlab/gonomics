package vcf

import "testing"

func TestToString(t *testing.T) {
	tests := []struct {
		record Vcf
		want   string
	}{
		{
			record: Vcf{
				Chr:    "chr1",
				Pos:    123456,
				Id:     "rs123456",
				Ref:    "A",
				Alt:    []string{"C", "G"},
				Qual:   99.9,
				Filter: "PASS",
				Info:   "DP=100",
			},
			want: "chr1\t123456\trs123456\tA\tC,G\t99.9\tPASS\tDP=100\n",
		},
		{
			record: Vcf{
				Chr:    "chr2",
				Pos:    234567,
				Id:     "rs234567",
				Ref:    "G",
				Alt:    []string{"A"},
				Qual:   50.5,
				Filter: "q10",
				Info:   "AF=0.5",
			},
			want: "chr2\t234567\trs234567\tG\tA\t50.5\tq10\tAF=0.5\n",
		},
	}

	for _, tt := range tests {
		t.Run(tt.record.Id, func(t *testing.T) {
			got := ToString(tt.record)
			if got != tt.want {
				t.Errorf("ToString() = %q, want %q", got, tt.want)
			}
		})
	}
}
