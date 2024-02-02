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
    contentAlpha := "chr9\t114159\t114165\tAhr::Arnt\t0\t+\t1\t0.2586407359495857\t0.7413592640504143\n"
    contentBeta := "chr9\t114159\t114165\tAhr::Arnt\t0\t+\t1\t0.2586407359495858\t0.7413592640504142\n"

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
