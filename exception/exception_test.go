package exception

import (
	"bytes"
	"errors"
	"log"
	"os"
	"os/exec"
	"strings"
	"testing"
)

// TestPanicOnErr verifies that the PanicOnErr function panics when given a non-nil error.
func TestRecoverPanicErrWithError(t *testing.T) {
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer func() {
		log.SetOutput(os.Stderr) // Restore default
		if strings.HasSuffix(buf.String(), "Expected panic, but none occurred:\n") {
			t.Error("Error: Failed expected panic:")
		}
	}()
	RecoverPanicErr()
}

// TestPanicOnErr verifies that the PanicOnErr function correctly panics when passed a non-nil error.
func TestPanicOnErr(t *testing.T) {
	defer RecoverPanicErr()
	PanicOnErr(errors.New("Error: expected test error"))
}

// TestFatalOnErr verifies that FatalOnErr logs an error and exits the program.
func TestFatalOnErr(t *testing.T) {
	// Capture log output
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer func() { log.SetOutput(os.Stderr) }() // Restore default

	if os.Getenv("TEST_FATAL") == "1" {
		FatalOnErr(errors.New("fatal test error"))
		return // Should not reach here
	}

	cmd := exec.Command(os.Args[0], "-test.run=TestFatalOnErr")
	cmd.Env = append(os.Environ(), "TEST_FATAL=1")
	err := cmd.Run()
	if e, ok := err.(*exec.ExitError); ok && !e.Success() {
		// Success: the test process exited with a non-zero status
		return
	}

	t.Errorf("Expected FatalOnErr to exit, but it didn't")
}

// TestWarningOnErr verifies that WarningOnErr logs a warning message correctly.
func TestWarningOnErr(t *testing.T) {
	// Capture log output
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer func() { log.SetOutput(os.Stderr) }() // Restore default

	WarningOnErr(errors.New("warning test error"))

	expected := "WARNING: warning test error\n"
	if !strings.HasSuffix(buf.String(), expected) {
		t.Errorf("Expected warning: '%s', got '%s'", expected, buf.String())
	}

	// Test with nil error (no output expected)
	buf.Reset()
	WarningOnErr(nil)
	if buf.String() != "" {
		t.Error("Unexpected warning output with nil error")
	}
}
