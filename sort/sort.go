package sort

import (
	"container/heap"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"os"
	"sort"
	"strconv"
	"sync"
)

const (
	maxTmpFilesAllowed = 1000
)

// TODO: I needed some way to determine which chunk the most recent giraf popped from
//  and came up with this. There may be a more elegant solution to this problem.
type priorityGiraf struct {
	data      *giraf.Giraf
	origin int
}

// Implements sort.Interface based on Breadth-First topological ordering of a SimpleGraph
type byTopologicalOrder []*priorityGiraf

func (g byTopologicalOrder) Len() int { return len(g) }

func (g byTopologicalOrder) Swap(i, j int) { g[i], g[j] = g[j], g[i] }

func (g byTopologicalOrder) Less(i, j int) bool {
	if g[i].data.Path.Nodes[0] < g[j].data.Path.Nodes[0] {
		return true
	} else if g[i].data.Path.Nodes[0] == g[j].data.Path.Nodes[0] {
		if g[i].data.Path.TStart < g[j].data.Path.TStart {
			return true
		}
	}
	return false
}

// Implements methods required for heap.Interface
func (g *byTopologicalOrder) Push(x interface{}) {
	giraf := x.(*priorityGiraf)
	*g = append(*g, giraf)
}

func (g *byTopologicalOrder) Pop() interface{} {
	oldQueue := *g
	n := len(oldQueue)
	giraf := oldQueue[n-1]
	oldQueue[n-1] = nil
	*g = oldQueue[:n-1]
	return giraf
}

// Sorts giraf by startposition and startnode given the node order defined in nodeIdSortOrder
// This was intended to be paired with a list of sorted nodes such as simpleGraph.GetSortOrder()
// so that the file result would be a topologically sorted giraf file, however this can be used with
// any given node sort order
func ExternalMergeSort(girafFile string, nodeIdSortOrder []uint32, linesPerChunk int, outFile string) {
	file := fileio.EasyOpen(girafFile)
	var done bool
	var chunkIDs []string
	var currChunkID string

	var currValues byTopologicalOrder
	sortOrderMap := sortOrderToMap(nodeIdSortOrder)

	for i := 0; !done; i++ {
		if i > maxTmpFilesAllowed {
			log.Fatalln("ERROR: Exceeded maximum number of tmp files, increase -chunkSize")
		}
		currValues = readChunk(file, linesPerChunk, &done)
		//TODO: is there a better way to ensure sorting by input nodeId order without renaming paths???
		switchIDs(currValues, sortOrderMap)
		sort.Sort(currValues)
		currChunkID = writeChunk(currValues, i)
		chunkIDs = append(chunkIDs, currChunkID)
	}
	file.Close()

	invertedSortOrder := invertMap(sortOrderMap)
	var wg sync.WaitGroup
	writeChan := make(chan *giraf.Giraf)

	go mergeChunks(writeChan, chunkIDs, invertedSortOrder)
	wg.Add(1)
	giraf.GirafChanToFile(outFile, writeChan, &wg)
	writeIdx(outFile, nodeIdSortOrder)
}

func sortOrderToMap(desiredOrder []uint32) map[uint32]uint32 {
	answer := make(map[uint32]uint32)
	var i uint32
	for i = 0; int(i) < len(desiredOrder); i++ {
		answer[i] = desiredOrder[i]
	}
	return answer
}

func invertMap(m map[uint32]uint32) map[uint32]uint32 {
	n := make(map[uint32]uint32)
	for key, value := range m {
		n[value] = key
	}
	return n
}

func switchIDs(g []*priorityGiraf, m map[uint32]uint32) {
	for i := 0; i < len(g); i++ {
		g[i].data.Path.Nodes[0] = m[g[i].data.Path.Nodes[0]]
	}
}

func readChunk(file *fileio.EasyReader, linesPerChunk int, done *bool) []*priorityGiraf {
	answer := make([]*priorityGiraf, linesPerChunk)
	var curr *giraf.Giraf
	for i := 0; i < linesPerChunk; i++ {
		curr, *done = giraf.NextGiraf(file)
		if *done {
			return answer[:i]
		}
		answer[i] = &priorityGiraf{curr, 0}
	}
	return answer
}

func writeChunk(g []*priorityGiraf, chunkNum int) string {
	chunkName := fmt.Sprintf("tmp_%d", chunkNum)
	file := fileio.MustCreate(chunkName)
	defer file.Close()
	for i := 0; i < len(g); i++ {
		giraf.WriteGirafToFileHandle(file, g[i].data)

	}
	return chunkName
}

func mergeChunks(outputChan chan<- *giraf.Giraf, chunkIDs []string, invertedSortOrder map[uint32]uint32) {
	var chunkReaders []*fileio.EasyReader = make([]*fileio.EasyReader, len(chunkIDs))
	priorityQueue := make(byTopologicalOrder, 0)
	// Gather the first element from each chunk to init heap
	for i := 0; i < len(chunkIDs); i++ {
		chunkReaders[i] = fileio.EasyOpen(chunkIDs[i])
		curr, done := giraf.NextGiraf(chunkReaders[i])
		if done {
			chunkReaders[i].Close()
			err := os.Remove(chunkReaders[i].File.Name())
			common.ExitIfError(err)
		} else {
			priorityQueue = append(priorityQueue, &priorityGiraf{curr, i})
		}
	}

	// Initialize min-heap
	heap.Init(&priorityQueue)

	var curr *priorityGiraf
	var nextGiraf *giraf.Giraf
	var done bool
	for priorityQueue.Len() > 0 {
		curr = heap.Pop(&priorityQueue).(*priorityGiraf)
		// Fix the path that we renamed to get the sort order
		curr.data.Path.Nodes[0] = invertedSortOrder[curr.data.Path.Nodes[0]]
		outputChan <- curr.data
		nextGiraf, done = giraf.NextGiraf(chunkReaders[curr.origin])
		if done {
			chunkReaders[curr.origin].Close()
			err := os.Remove(chunkReaders[curr.origin].File.Name())
			common.ExitIfError(err)
		} else {
			heap.Push(&priorityQueue, &priorityGiraf{nextGiraf, curr.origin})
		}
	}
	close(outputChan)
}

//TODO: should really come up with a better index format than this
func writeIdx(filename string, sortOrder []uint32) {
	idxName := filename + ".idx"
	file := fileio.MustCreate(idxName)
	defer file.Close()

	for i := 0; i < len(sortOrder); i++ {
		_, err := fmt.Fprintln(file, sortOrder[i])
		common.ExitIfError(err)
	}
}

func ReadIdx(filename string) []uint32 {
	answer := make([]uint32, 0)
	file := fileio.EasyOpen(filename)
	defer file.Close()

	for line, done := fileio.EasyNextLine(file); !done; line, done = fileio.EasyNextLine(file) {
		val, err := strconv.Atoi(line)
		common.ExitIfError(err)
		answer = append(answer, uint32(val))
	}
	return answer
}
