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

type priorityGiraf struct {
	data     *giraf.Giraf
	origin   int
	sortPath []uint32
}

// Implements sort.Interface based on Breadth-First topological ordering of a SimpleGraph
type byTopologicalOrder []*priorityGiraf

func (g byTopologicalOrder) Len() int { return len(g) }

func (g byTopologicalOrder) Swap(i, j int) { g[i], g[j] = g[j], g[i] }

func (g byTopologicalOrder) Less(i, j int) bool {
	// determine minPathLen so that we dont overshoot the index on one of the inputs
	var minPathLen int
	if len(g[i].sortPath) < len(g[j].sortPath) {
		minPathLen = len(g[i].sortPath)
	} else {
		minPathLen = len(g[j].sortPath)
	}

	// check first node and see if there is a sortable difference
	if g[i].sortPath[0] < g[j].sortPath[0] {
		return true
	} else if g[i].sortPath[0] > g[j].sortPath[0] {
		return false
		// if the nodes are the same, check the start position
	} else {
		if g[i].data.Path.TStart < g[j].data.Path.TStart {
			return true
		} else if g[i].data.Path.TStart > g[j].data.Path.TStart {
			return false
			// if the start pos is the same go through the paths until there is a sortable difference
		} else {
			for k := 1; k < minPathLen-1; k++ {
				if g[i].sortPath[k] < g[j].sortPath[k] {
					return true
				} else if g[i].sortPath[k] > g[j].sortPath[k] {
					return false
				}
			}
		}
	}

	// if there is no sortable difference, the shorter path is less
	if len(g[i].sortPath) < len(g[j].sortPath) {
		return true
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
	//oldQueue[n-1] = nil // potential solution to memory leak?
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
		currValues = readChunk(file, linesPerChunk, &done, sortOrderMap)
		// Check to see if currValues contains any data. If not, continue.
		if len(currValues) == 0 {
			continue
		}
		sort.Sort(currValues)
		currChunkID = writeChunk(currValues, i)
		chunkIDs = append(chunkIDs, currChunkID)
	}
	file.Close()

	//invertedSortOrder := invertMap(sortOrderMap)
	var wg sync.WaitGroup
	writeChan := make(chan *giraf.Giraf)

	go mergeChunks(writeChan, chunkIDs, sortOrderMap)
	wg.Add(1)
	giraf.GirafChanToFile(outFile, writeChan, &wg)
	writeIdx(outFile, nodeIdSortOrder)
}

func sortOrderToMap(desiredOrder []uint32) map[uint32]uint32 {
	answer := make(map[uint32]uint32)
	var i uint32
	for i = 0; int(i) < len(desiredOrder); i++ {
		answer[desiredOrder[i]] = i
	}
	return answer
}

func getSortPath(path *giraf.Path, sortOrderMap map[uint32]uint32) []uint32 {
	answer := make([]uint32, len(path.Nodes))
	for i := 0; i < len(path.Nodes); i++ {
		answer[i] = sortOrderMap[path.Nodes[i]]
	}
	return answer
}

func readChunk(file *fileio.EasyReader, linesPerChunk int, done *bool, sortOrderMap map[uint32]uint32) []*priorityGiraf {
	answer := make([]*priorityGiraf, 0, linesPerChunk)
	var curr *giraf.Giraf
	for i := 0; i < linesPerChunk; i++ {
		curr, *done = giraf.NextGiraf(file)
		if *done {
			return answer
		}
		answer = append(answer, &priorityGiraf{curr, 0, getSortPath(curr.Path, sortOrderMap)})
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

func mergeChunks(outputChan chan<- *giraf.Giraf, chunkIDs []string, sortOrderMap map[uint32]uint32) {
	var chunkReaders []*fileio.EasyReader = make([]*fileio.EasyReader, len(chunkIDs))
	priorityQueue := make(byTopologicalOrder, 0, len(chunkIDs))
	var sortPath []uint32
	// Gather the first element from each chunk to init heap
	for i := 0; i < len(chunkIDs); i++ {
		chunkReaders[i] = fileio.EasyOpen(chunkIDs[i])
		curr, done := giraf.NextGiraf(chunkReaders[i])
		if done {
			chunkReaders[i].Close()
			err := os.Remove(chunkReaders[i].File.Name())
			common.ExitIfError(err)
		} else {
			// I figure that it would be just as quick to rederive the sortPath
			// compared to writing it to the tmp file and decoding it here
			sortPath = getSortPath(curr.Path, sortOrderMap)
			priorityQueue = append(priorityQueue, &priorityGiraf{curr, i, sortPath})
		}
	}

	// Initialize min-heap
	heap.Init(&priorityQueue)

	var curr *priorityGiraf
	var nextGiraf *giraf.Giraf
	var done bool
	for priorityQueue.Len() > 0 {
		curr = heap.Pop(&priorityQueue).(*priorityGiraf)
		outputChan <- curr.data
		nextGiraf, done = giraf.NextGiraf(chunkReaders[curr.origin])
		if done {
			chunkReaders[curr.origin].Close()
			err := os.Remove(chunkReaders[curr.origin].File.Name())
			common.ExitIfError(err)
		} else {
			heap.Push(&priorityQueue, &priorityGiraf{nextGiraf, curr.origin, getSortPath(nextGiraf.Path, sortOrderMap)})
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
