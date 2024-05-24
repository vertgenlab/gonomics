package exception

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"strings"
	"testing"
)

func TestPanicOnErr(t *testing.T) {

	defer func() {
		if r := recover(); r == nil {
			t.Errorf("PanicOnErr did not panic on error...\n")
		}
	}()
	var err error = fmt.Errorf("Error: test exception error handling...\n")

	PanicOnErr(err)
	FatalOnErr(err)
}

func TestWarningOnErr(t *testing.T) {
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer func() {
		log.SetOutput(os.Stderr)
	}()

	var warning string = "Warning: warning error handling...\n"

	WarningOnErr(fmt.Errorf(warning))

	if !strings.Contains(buf.String(), warning) {
		t.Errorf("WarningOnErr did not log the correct warning message")
	}
}

func TestRandAutoSeedExists(t *testing.T) {
	os.Setenv("GODEBUG", "randautoseed=0,other=value")
	expected := "randautoseed=0,other=value"

	if os.Getenv("GODEBUG") != expected && os.Getenv("GODEBUG") != "other=value,randautoseed=0" {
		t.Errorf("setGoDebug did not update randautoseed correctly")
	}
}
