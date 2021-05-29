package main

import (
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"runtime"
	"sort"
	"strings"
)

var Reset  = "\033[0m"
var Red    = "\033[31m"
var Green  = "\033[32m"
var Yellow = "\033[33m"
var Blue   = "\033[34m"
var Purple = "\033[35m"
var Cyan   = "\033[36m"
var Gray   = "\033[37m"
var White  = "\033[97m"

func init() {
	if runtime.GOOS == "windows" {
		Reset  = ""
		Red    = ""
		Green  = ""
		Yellow = ""
		Blue   = ""
		Purple = ""
		Cyan   = ""
		Gray   = ""
		White  = ""
	}
}

func usage() {
	cache := getCache()
	fmt.Print(
		"gonomics - A collection of builtin tools that use the gonomics core library.\n\n" +
			"Usage: gonomics <command> [options]\n\n")
	printCmdList(cache)
}

// getCache returns a cache file ('go bin path'/.cmdcache) listing the cmd usage statements.
// Retrieves the file if it exists, else builds a new cache.
func getCache() (cache *os.File){
	binPath, _ := getBin()

	cmdStat, _ := os.Stat(os.Args[0]) // stat binary for this file
	cacheStat, cacheStatErr := os.Stat(binPath + "/.cmdcache") // stat cache file

	// Build a cache if cache does not exist OR if the cache is older than the gonomics cmd binary
	if os.IsNotExist(cacheStatErr) || cmdStat.ModTime().After(cacheStat.ModTime()) {
		buildCmdCache(binPath)
	}
	file, err := os.Open(binPath + "/.cmdcache")
	if err != nil {
		log.Panic(err)
	}
	return file
}

// printCmdList prints the cache file
func printCmdList(cache *os.File) {
	toPrint, err := io.ReadAll(cache)
	cache.Close()
	if err != nil {
		log.Panic(err)
	}
	fmt.Print(string(toPrint))
}

func buildCmdCache(binPath string) {
	cacheWriter, err := os.Create(binPath + "/.cmdcache")
	defer cacheWriter.Close()
	if err != nil {
		log.Panic(err)
	}

	cmdMap := getGonomicsCmds()

	cmds := make([]string, 0, len(cmdMap))
	for key := range cmdMap {
		cmds = append(cmds, key)
	}
	sort.Slice(cmds, func(i, j int) bool { return cmds[i] < cmds[j] })

	var endFirstLineIdx int
	var rawOutput []byte

	_, err = fmt.Fprintln(cacheWriter, Red+"Commands:"+Reset)
	if err != nil {
		log.Panic(err)
	}

	for _, cmdName := range cmds {
		if cmdName == "gonomics" { // avoid recursive call of the gonomics cmd
			continue
		}

		cmd := exec.Command(binPath + "/" + cmdName)
		rawOutput, _ = cmd.Output()
		endFirstLineIdx = strings.Index(string(rawOutput), "\n")

		switch {
		case endFirstLineIdx > 0: // cmd starts with summary line
			_, err = fmt.Fprintf(cacheWriter, "     %s\n", string(rawOutput[:endFirstLineIdx]))

		case len(rawOutput) == 0: // cmd has no usage statement
			_, err = fmt.Fprintf(cacheWriter, "     %s\n", cmdName)

		case endFirstLineIdx == -1 && len(rawOutput) > 0: // has output, but no newline
			_, err = fmt.Fprintf(cacheWriter, "     %s\n", rawOutput)

		default: // if all else fails, print cmd name
			_, err = fmt.Fprintf(cacheWriter, "     %s\n", cmdName)
		}

		if err != nil {
			log.Panic(err)
		}
	}
}
