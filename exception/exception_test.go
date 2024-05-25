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
	var err error = fmt.Errorf("Error: panic exception error handling...\n")
	PanicOnErr(err)
}

func TestFatalOnErr(t *testing.T) {
	// Test with no error
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer func() { log.SetOutput(os.Stderr) }()

	FatalOnErr(nil)
	if buf.Len() > 0 {
		t.Errorf("Error: FatalOnErr logged output when no error was provided")
	}

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
		t.Errorf("Error: WarningOnErr did not log the correct warning message")
	}
}
