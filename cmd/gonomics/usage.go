//# GONOMICS: testse

// Package main does stuff
package main

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"runtime"
	"sort"
	"strings"
	"text/tabwriter"
)

var Reset = "\033[0m"
var Red = "\033[31m"
var Green = "\033[32m"
var Yellow = "\033[33m"
var Blue = "\033[34m"
var Purple = "\033[35m"
var Cyan = "\033[36m"
var Gray = "\033[37m"
var White = "\033[97m"

const hline = "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

// init prevents color characters if running Windows.
// Windows does not allow colored terminal output by
// default so it just adds a bunch of characters to
// the print statements.
func init() {

	if runtime.GOOS == "windows" {
		Reset = ""
		Red = ""
		Green = ""
		Yellow = ""
		Blue = ""
		Purple = ""
		Cyan = ""
		Gray = ""
		White = ""
	}
}

type CmdInfo struct {
	Name  string
	Usage string
	Group string
}

func usage() {
	cache := getCache()

	fmt.Print(
		Yellow + "gonomics - A collection of tools that use the gonomics core library.\n\n" +
			"Usage: gonomics <command> [options]\n\n" + Reset)
	fmt.Println(Red + "Commands:" + Reset)
	fmt.Print(hline)
	printCmdList(cache)
	fmt.Println(hline)
}

// getCache returns a cache file ('go bin path'/.cmdcache) listing the cmd usage statements.
// Retrieves the file if it exists, else builds a new cache.
func getCache() (cache *os.File) {
	binPath, _ := getBin()
	cmdStat, _ := os.Stat(binPath + "/gonomics")               // stat binary for this file
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
	toPrint, err := ioutil.ReadAll(cache)
	cache.Close()
	if err != nil {
		log.Panic(err)
	}
	fmt.Print(string(toPrint))
}

// buildCmdCache creates a cache file listing the cmd usage statements.
// TODO find a way to store cache in gonomics binary file
func buildCmdCache(binPath string) {
	fmt.Println("...Building Cache...")
	cmdMap := getGonomicsCmds()

	cmds := make([]string, 0, len(cmdMap))
	for key := range cmdMap {
		cmds = append(cmds, key)
	}
	sort.Slice(cmds, func(i, j int) bool { return cmds[i] < cmds[j] })

	var cmdSourcePath string
	var cachedSrcPath string = getCachedSrcDir(binPath + "/.cmdcache")
	switch {
	case cachedSrcPath != "":
		cmdSourcePath = cachedSrcPath
	case os.Getenv("GOPATH") != "":
		cmdSourcePath = os.Getenv("GOPATH") + "/src/github.com/vertgenlab/gonomics/cmd/"
	default:
		cmdSourcePath = os.Getenv("HOME") + "/src/github.com/vertgenlab/gonomics/cmd/"
	}

	var cmdGroup, cmdUsage string

	groupMap := make(map[string][]CmdInfo) // key: group, value: all cmds in group

	for _, cmdName := range cmds {
		if cmdName == "gonomics" { // avoid recursive call of the gonomics cmd
			continue
		}

		// Open source file and parse gonomics cmd tags
		cmdGroup, cmdUsage = parseTagsFromSource(cmdSourcePath + cmdName + "/" + cmdName + ".go")

		if cmdGroup == "" { // TODO grouping map
			cmdGroup = "Uncategorized"
		}

		if cmdUsage == "" {
			cmdUsage = getUsageFromRunningCmd(binPath + "/" + cmdName)
		}

		info := CmdInfo{Name: cmdName, Usage: cmdUsage, Group: cmdGroup}

		groupMap[info.Group] = append(groupMap[info.Group], info)
	}
	writeCache(binPath, groupMap, cmdSourcePath)
}

// writeCache writes a groupMap to file
// TODO change so cache is added to gonomics cmd binary
func writeCache(binPath string, groupMap map[string][]CmdInfo, srcPath string) {
	cacheWriter, err := os.Create(binPath + "/.cmdcache")
	defer cacheWriter.Close()
	if err != nil {
		log.Panic(err)
	}

	// write source directory line to cache
	_, err = fmt.Fprintf(cacheWriter, "##SourceDir:%s\n", srcPath)
	if err != nil {
		log.Panic(err)
	}

	// initialize tabwriter
	w := new(tabwriter.Writer)

	// minwidth, tabwidth, padding, padchar, flags
	w.Init(cacheWriter, 4, 8, 0, ' ', 0)

	defer w.Flush()

	var allGroups [][]CmdInfo
	for _, cmds := range groupMap {
		allGroups = append(allGroups, cmds)
	}

	sort.Slice(allGroups, func(i, j int) bool {
		if allGroups[i][0].Group == "Uncategorized" { // want this at the bottom
			return false
		}
		if allGroups[j][0].Group == "Uncategorized" { // want this at the bottom
			return true
		}
		return allGroups[i][0].Group < allGroups[j][0].Group
	})

	for _, cmds := range allGroups {
		_, err = fmt.Fprintf(w, "\t\n%s\t\t\n", Red+cmds[0].Group+Reset)
		if err != nil {
			log.Panic(err)
		}

		for _, cmd := range cmds {
			_, err = fmt.Fprintf(w, "     %s\t%s\n", Green+cmd.Name+Reset, Cyan+cmd.Usage+Reset)
			if err != nil {
				log.Panic(err)
			}
		}
	}
}

// parseTagsFromSource parses gonomics cmd tags from source files.
func parseTagsFromSource(filepath string) (group string, usage string) {
	tagLines := getHeaderCommentLines(filepath)

	for _, line := range tagLines {
		switch {
		case strings.Contains(line, "Command Group:"):
			group = strings.FieldsFunc(line, func(c rune) bool { return c == '"' })[1]
		case strings.Contains(line, "Command Usage:"):
			usage = strings.FieldsFunc(line, func(c rune) bool { return c == '"' })[1]
		}
	}

	return
}

// getUsageFromRunningCmd retrieves the usage statement by parsing the
// output of the command.
func getUsageFromRunningCmd(cmdPath string) (usage string) {
	var rawOutput []byte
	cmd := exec.Command(cmdPath)
	rawOutput, _ = cmd.Output()
	endFirstLineIdx := strings.Index(string(rawOutput), "\n")
	hypenIdx := strings.Index(string(rawOutput), "-")

	if hypenIdx > endFirstLineIdx {
		hypenIdx = 0
	}

	switch {
	case endFirstLineIdx > 0: // cmd starts with summary line
		return strings.TrimSpace(string(rawOutput[hypenIdx+1 : endFirstLineIdx]))

	case len(rawOutput) == 0: // cmd has no usage statement
		return ""

	case endFirstLineIdx == -1 && len(rawOutput) > 0: // has output, but no newline
		return strings.TrimSpace(string(rawOutput[hypenIdx+1:]))

	default: // if all else fails, print cmd name
		return ""
	}
}

// getHeaderCommentLines retrieves all header lines in a file that begin
// with the string "//" (golang comment string).
func getHeaderCommentLines(filepath string) []string {
	var answer []string
	var line string
	var done bool
	file := fileio.EasyOpen(filepath)

	for line, done = fileio.EasyNextLine(file); !done && strings.HasPrefix(line, "//"); line, done = fileio.EasyNextLine(file) {
		answer = append(answer, line)
	}
	err := file.Close()
	if err != nil {
		log.Panic(err)
	}
	return answer
}

// getCachedSrcDir retrieves the source code directory stored in the first line of the cache file
func getCachedSrcDir(cacheFile string) string {
	file, err := os.Open(cacheFile)
	defer file.Close()
	if err == os.ErrNotExist {
		return ""
	} else if err != nil {
		log.Panic(err)
	}

	s := bufio.NewScanner(file)
	s.Scan()
	line := s.Text()
	if strings.HasPrefix(line, "##SourceDir:") {
		return strings.TrimPrefix(line, "##SourceDir:")
	}
	exception.PanicOnErr(err)
	return ""
}

// buildFromPath is run when setpath is called and write a new cache file with specified gonomicsPath
func buildFromPath(gonomicsPath, binPath string) {
	gonomicsPath = strings.TrimRight(gonomicsPath, "/") + "/cmd/"
	file, err := os.Create(binPath + "/.cmdcache")
	defer file.Close()
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(file, "##SourceDir:%s\n", gonomicsPath)
	if err != nil {
		log.Panic(err)
	}
}
