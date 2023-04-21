package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestGonomics(t *testing.T) {
	// check to see if binary is already installed, if so buildCmdCache, else return

	// Find binary location. Preference is gonomics exec path > GOBIN > GOPATH > Default go install location
	execPath, err := os.Executable() // path to gonomics executable
	exception.PanicOnErr(err)
	execPath, err = filepath.EvalSymlinks(execPath) // follow any symlinks to get true executable location

	// options to choose from
	execPath = strings.TrimSuffix(execPath, ".test") // avoid errors when testing function
	execPath = strings.TrimSuffix(execPath, "/gonomics")
	gobin := os.Getenv("GOBIN")
	gopath := os.Getenv("GOPATH") + "/bin"     // default to $GOPATH/bin
	godefault := os.Getenv("HOME") + "/go/bin" // default install directory

	switch {
	case tryPathFile(execPath):
	case gobin != "" && tryPathFile(gobin):
	case gopath != "" && tryPathFile(gopath):
	case godefault != "" && tryPathFile(godefault):
	default:
		t.SkipNow()
	}
	buildCmdCache("") // tests everything is formatted correctly and updates cache file
}
