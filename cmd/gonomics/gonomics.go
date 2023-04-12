package main

import (
	_ "embed"
	"flag"
	"github.com/vertgenlab/gonomics/exception"
	"io/fs"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

//go:embed command_cache.txt
var cache string

// getBin returns the path to binary executables and a map of existing executables
func getBin() (path string, binExists map[string]bool) {
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
		path = execPath

	case gobin != "" && tryPathFile(gobin):
		path = gobin

	case gopath != "" && tryPathFile(gopath):
		path = gopath

	case godefault != "" && tryPathFile(godefault):
		path = godefault

	default:
		log.Fatalf("ERROR: could not find executables for gonomics in any of the following paths:"+
			"\n'%s'\n'%s'\n'%s'\n'%s'\n", execPath, gobin, gopath, godefault)
	}

	var dir []fs.FileInfo
	dir, err = ioutil.ReadDir(path)
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

// getSrc returns the path to the source code files for commands
func getSrc() string {
	var expectedPath string
	var cachedSrcPath string = getCachedSrcDir()
	gopath := os.Getenv("GOPATH") + "/src/github.com/vertgenlab/gonomics/cmd/"
	godefault := os.Getenv("HOME") + "/go/src/github.com/vertgenlab/gonomics/cmd/"
	currwd, err := os.Getwd()
	exception.PanicOnErr(err)

	switch {
	case cachedSrcPath != "" && tryPathDir(cachedSrcPath):
		expectedPath = cachedSrcPath

	case os.Getenv("GOPATH") != "" && tryPathDir(gopath):
		expectedPath = gopath

	case os.Getenv("HOME") != "" && tryPathDir(godefault):
		expectedPath = godefault

	case tryPathDir(currwd):
		expectedPath = currwd

	default:
		log.Fatalf("ERROR: could not find gonomics cmd folder in any of the following locations\n%s\n%s\n%s\n"+
			"Please use the '-setpath' flag followed by the path to the gonomics directory\n"+
			"Subsequent calls of the gonomics command will not require the '-setpath' flag\n", cachedSrcPath, gopath, godefault)
	}
	return expectedPath
}

// tryPathFile tests a path for the presence of intervalOverlap (arbitrary chosen command)
// and returns true if the intervalOverlap executable is found in the input path.
func tryPathFile(path string) bool {
	var err error
	var dir []fs.FileInfo
	_, err = os.Stat(path)
	if err != nil {
		return false
	}
	dir, err = ioutil.ReadDir(path)
	if err != nil {
		log.Printf("ERROR: could not access %s\n", path)
		log.Fatal(err)
	}
	for _, file := range dir {
		if file.Name() == "intervalOverlap" && !file.IsDir() {
			return true
		}
	}
	return false
}

// tryPathDir tests a path for the presence of intervalOverlap source code directory
// and returns true if the intervalOverlap directory is found in the input path.
func tryPathDir(path string) bool {
	var err error
	var dir []fs.FileInfo
	_, err = os.Stat(path)
	if err != nil {
		return false
	}
	dir, err = ioutil.ReadDir(path)
	if err != nil {
		log.Printf("ERROR: could not access %s\n", path)
		log.Fatal(err)
	}
	for _, file := range dir {
		if file.Name() == "intervalOverlap" && file.IsDir() {
			return true
		}
	}
	return false
}

// getGonomicsCmds parses the gonomics source code to return a set of cmd names
// path input pointing to gonomics directory overrides autodetection
func getGonomicsCmds(altSrcPath string) map[string]bool {
	if cache != "" {
		ans := getCmdsFromCache()
		if ans["gonomics"] { // check that the directory for this cmd was found
			return ans
		}
	} // else try to parse from src files

	var srcPath string
	if altSrcPath != "" {
		srcPath = altSrcPath
	} else {
		srcPath = getSrc()
	}
	cmds, err := ioutil.ReadDir(srcPath)
	if err != nil {
		log.Fatalf("ERROR: could not find gonomics cmd folder in expected path: %s\n"+
			"Please use the '-setpath' flag followed by the path to the gonomics directory\n"+
			"Subsequent calls of the gonomics command will not require the '-setpath' flag\n", srcPath)
	}

	funcNames := make(map[string]bool)

	for _, x := range cmds {
		funcNames[x.Name()] = true
	}

	if !funcNames["gonomics"] { // check that the directory for this cmd was found
		log.Fatalf("ERROR: gonomics cmd folder was not found in %s\n", srcPath)
	}

	return funcNames
}

// getCmdsFromCache parses the cached file in and returns the cmd map.
func getCmdsFromCache() map[string]bool {
	funcNames := make(map[string]bool)
	lines := strings.Split(cache, "\n")
	for i := range lines {
		if strings.HasPrefix(lines[i], "##CmdName:") {
			funcNames[strings.TrimPrefix(lines[i], "##CmdName:")] = true
		}
	}
	return funcNames
}

func main() {
	gonomicsPath := flag.String("setpath", "", "path to gonomics directory")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *gonomicsPath != "" {
		buildCmdCache(*gonomicsPath)
		return
	}

	if len(flag.Args()) == 0 {
		flag.Usage()
		return
	}

	cmdCalled := flag.Arg(0) // command called by the user (e.g. 'gonomics faFormat')

	binPath, binExists := getBin()
	isGonomicsCmd := getGonomicsCmds(*gonomicsPath)

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
