//# GONOMICS: testse

// Package main does stuff
package main

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"runtime"
	"sort"
	"strings"
	"text/tabwriter"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"golang.org/x/exp/slices"
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
	getCache()

	fmt.Print(
		Yellow + "gonomics - A collection of tools that use the gonomics core library.\n\n" +
			"Usage: gonomics <command> [options]\n\n" + Reset)
	fmt.Println(Red + "Commands:" + Reset)
	fmt.Print(hline)
	printCmdList()
	fmt.Println(hline)
}

// getCache makes sure the cache is present, and builds a new cache if it is not.
func getCache() {
	if cache != "" {
		return // cache is already present
	}

	buildCmdCache("")
}

// printCmdList prints the cache file.
func printCmdList() {
	lines := strings.Split(cache, "\n")
	var lastIdx int
	for i := range lines {
		if strings.HasPrefix(lines[i], "##") {
			lastIdx = i
		}
	}
	fmt.Print(strings.Join(lines[lastIdx+1:], "\n"))
}

// buildCmdCache builds the cache listing the cmd usage statements.
func buildCmdCache(altSrcPath string) {
	log.Println("...Building Cache...")

	// get cmd names and manipulate to sorted slice
	cmdMap := getGonomicsCmds(altSrcPath)
	cmds := make([]string, 0, len(cmdMap))
	for key := range cmdMap {
		cmds = append(cmds, key)
	}
	sort.Slice(cmds, func(i, j int) bool { return cmds[i] < cmds[j] })
	// get location of binary and source files
	binPath, _ := getBin()

	var srcPath string
	if altSrcPath != "" {
		srcPath = altSrcPath
	} else {
		srcPath = getSrc()
	}

	// sort cmds into command groups
	var cmdGroup, cmdUsage string
	groupMap := make(map[string][]CmdInfo) // key: group, value: all cmds in group

	for _, cmdName := range cmds {
		if cmdName == "gonomics" { // avoid recursive call of the gonomics cmd
			continue
		}

		// Open source file and parse gonomics cmd tags
		cmdGroup, cmdUsage = parseTagsFromSource(srcPath + cmdName + "/" + cmdName + ".go")
		if cmdGroup == "" { // TODO grouping map
			cmdGroup = "Uncategorized"
		}
		if cmdUsage == "" {
			cmdUsage = getUsageFromRunningCmd(binPath + "/" + cmdName)
		}

		info := CmdInfo{Name: cmdName, Usage: cmdUsage, Group: cmdGroup}

		groupMap[info.Group] = append(groupMap[info.Group], info)
	}

	writeCache(groupMap, srcPath)
}

// writeCache writes a groupMap to file.
func writeCache(groupMap map[string][]CmdInfo, srcPath string) {
	outfile := srcPath + "/gonomics/command_cache.txt"
	cacheWriter, err := os.Create(outfile)
	if err != nil {
		log.Panic(err)
	}

	// write source directory line to cache
	_, err = fmt.Fprintf(cacheWriter, "##SourceDir:%s\n", srcPath)
	if err != nil {
		log.Panic(err)
	}

	// write command names
	var cmdNames []string
	for _, group := range groupMap {
		for i := range group {
			cmdNames = append(cmdNames, fmt.Sprintf("##CmdName:%s\n", group[i].Name))
		}
	}
	slices.Sort(cmdNames) // so file is not updated on PRs with no cmd changes
	for i := range cmdNames {
		_, err = fmt.Fprintf(cacheWriter, cmdNames[i])
		exception.PanicOnErr(err)
	}

	// initialize tabwriter
	w := new(tabwriter.Writer)

	// minwidth, tabwidth, padding, padchar, flags
	w.Init(cacheWriter, 4, 8, 0, ' ', 0)

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

	err = w.Flush()
	exception.PanicOnErr(err)
	err = cacheWriter.Close()
	exception.PanicOnErr(err)

	cache = strings.Join(fileio.Read(outfile), "\n") + "\n"
}

// parseTagsFromSource parses gonomics cmd tags from source files.
func parseTagsFromSource(filepath string) (group string, usage string) {
	tagLines := getHeaderCommentLines(filepath)

	// ATTENTION //
	// If your PR is failing tests due to the gonomics command tests, please ensure that if you declare a
	// command group for your cmd, you do so in the following format (ignore first // on each line):
	// // Command Group: "Deep Learning"
	//
	// package main
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

// getCachedSrcDir retrieves the source code directory stored in the first line of the cache.
func getCachedSrcDir() string {
	lines := strings.Split(cache, "\n")
	for i := range lines {
		if strings.HasPrefix(lines[i], "##SourceDir:") {
			return strings.TrimPrefix(lines[i], "##SourceDir:")
		}
	}
	return ""
}
