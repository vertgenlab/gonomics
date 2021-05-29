package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/exception"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
)

func getBin() (path string, binExists map[string]bool) {
	switch { // Find binary location. Preference is GOBIN > GOPATH > Default go install location
	case os.Getenv("GOBIN") != "":
		path = os.Getenv("GOBIN")

	case os.Getenv("GOPATH") != "":
		path = os.Getenv("GOPATH") + "/bin" // default to $GOPATH/bin

	default:
		path = os.Getenv("HOME") + "/go/bin" // default install directory
	}

	dir, err := ioutil.ReadDir(path)
	if err != nil {
		log.Printf("ERROR: could not access %s\n", path)
		log.Fatal(err)
	}

	binExists = make(map[string]bool)

	for _, file := range dir {
		binExists[file.Name()] = true
	}

	if !binExists["gonomics"] { // check if gonomics is installed by checking for this binary
		log.Fatalf("ERROR: gonomics binary were not found in %s\n "+
			"Run 'go install ./...' in gonomics directory to compile binary", path)
	}

	return path, binExists
}

func getGonomicsCmds() map[string]bool {
	expectedPath := os.Getenv("GOPATH") + "/src/github.com/vertgenlab/gonomics/cmd"
	cmds, err := ioutil.ReadDir(expectedPath)
	if err != nil {
		log.Printf("ERROR: could not find gonomics cmd folder in expected path: %s\n", expectedPath)
		log.Fatal(err)
	}

	funcNames := make(map[string]bool)

	for _, x := range cmds {
		funcNames[x.Name()] = true
	}

	if !funcNames["gonomics"] { // check that the directory for this cmd was found
		log.Fatalf("ERROR: gonomics cmd folder was not found in %s\n", expectedPath)
	}

	return funcNames
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) == 0 {
		flag.Usage()
		return
	}

	cmdCalled := flag.Arg(0) // command called by the user (e.g. 'gonomics faFormat')

	binPath, binExists := getBin()
	isGonomicsCmd := getGonomicsCmds()

	switch { // Error checks. verbose for clarity
	case isGonomicsCmd[cmdCalled] && binExists[cmdCalled]:
		// proceed in main function

	case isGonomicsCmd[cmdCalled] && !binExists[cmdCalled]:
		log.Fatalf("ERROR: binary for %s was not found. "+
			"Please run 'go install ./...' in the gonomics directory.", cmdCalled)

	case !isGonomicsCmd[cmdCalled] && binExists[cmdCalled]:
		log.Fatalf("ERROR: %s is not a gonomics command", cmdCalled)

	default:
		log.Fatalf("ERROR: %s does not exist, and is not a gonomics command", cmdCalled)
	}

	cmd := exec.Command(binPath+"/"+cmdCalled, flag.Args()[1:]...)

	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	exception.FatalOnErr(err)
}
