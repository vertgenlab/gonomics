package motif

import (
	"io/ioutil"
	"os"
	"testing"
)

func createTempFileWithContent(t *testing.T, content string) (filename string) {
	t.Helper()
	tmpFile, err := ioutil.TempFile("", "motif_test_*.tsv")
	if err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	defer tmpFile.Close()

	_, err = tmpFile.WriteString(content)
	if err != nil {
		t.Fatalf("Failed to write to temp file: %v", err)
	}
	return tmpFile.Name()
}

func TestApproxEquals(t *testing.T) {
	epsilon := 1e-12
	contentAlpha := "chr1\t1\t100\tGeneX\t0\t+\t1\t0.2586407359495857\t0.7413592640504143\n"
	contentBeta := "chr1\t1\t100\tGeneX\t0\t+\t1\t0.2586407359495858\t0.7413592640504142\n"

	alphaPath := createTempFileWithContent(t, contentAlpha)
	betaPath := createTempFileWithContent(t, contentBeta)
	defer os.Remove(alphaPath)
	defer os.Remove(betaPath)

	got := ApproxEquals(alphaPath, betaPath, epsilon)
	want := true
	if got != want {
		t.Errorf("ApproxEquals() = %v; want %v", got, want)
	}
}

func TestApproxEqualsExtended(t *testing.T) {
	epsilon := 1e-12
	tests := []struct {
		name         string
		contentAlpha string
		contentBeta  string
		want         bool
	}{
		{
			name:         "DifferingNumberOfLines",
			contentAlpha: "chrX\t114159\t114165\tGeneX\t0\t+\t1\t0.2586407359495857\t0.7413592640504143\n",
			contentBeta: "chrX\t114159\t114165\tGeneX\t0\t+\t1\t0.2586407359495857\t0.7413592640504143\n" +
				"chr9\t114050\t114056\tZNF354C\t0\t-\t0\t0.8165285748498434\t0.8165285748498434\n",
			want: false,
		},
		{
			name:         "DifferingNumberOfFields",
			contentAlpha: "chrX\t114159\t114165\tGeneX\t0\t+\t1\t0.2586407359495857\n",
			contentBeta:  "chr9\t114159\t114165\tGeneX\t0\t+\t1\n",
			want:         false,
		},
		{
			name:         "IndexOutOfRange",
			contentAlpha: "chrX\t114159\t114165\tGeneX\t0\t+\t1\n",
			contentBeta:  "chrX\t114159\t114165\tGeneX\t0\t+\t1\n",
			want:         false,
		},
		{
			name:         "NotEqual",
			contentAlpha: "chr1\t0\t100\tGeneX\t0\t+\t1\t0.2586407359495857\t0.7413592640504143\n",
			contentBeta:  "chr1\t0\t100\tGeneX\t0\t+\t1\t0.3586407359495857\t0.8413592640504143\n",
			want:         false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			alphaPath := createTempFileWithContent(t, tt.contentAlpha)
			betaPath := createTempFileWithContent(t, tt.contentBeta)
			defer os.Remove(alphaPath)
			defer os.Remove(betaPath)

			got := ApproxEquals(alphaPath, betaPath, epsilon)
			if got != tt.want {
				t.Errorf("ApproxEquals() for %s = %v; want %v", tt.name, got, tt.want)
			}
		})
	}
}
