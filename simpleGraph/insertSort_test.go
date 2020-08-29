package simpleGraph

/*
import(
	"sort"
	"testing"
	"log"
)

var seeds []int = []int{4, 1, 3, 2, 16, 9, 10, 14, 8, 7}

func TestInsertSort(t *testing.T) {

	log.Printf("%v\n", seeds)


	heapSorted := heapSort(seeds)
	log.Printf("%v\n", heapSorted)
}

func insertSort(data []int, el int) []int {
	index := sort.Search(len(data), func(i int) bool { return data[i] <= el })
	data = append(data, 0)
	copy(data[index+1:], data[index:])
	data[index] = el
	return data
}

func insertSortAll(answer []*pathSeed, newData []*pathSeed) []*pathSeed {
	finalNum := len(answer)
	answer = append(answer, make([]*pathSeed, len(newData))...)
	var index int
	for i := 0; i < len(newData); i++ {
		index = sort.Search(finalNum, func(j int) bool { return answer[j].Len <= newData[i].Len })
		if index < finalNum {
			copy(answer[index+1:], answer[index:])
		}
		answer[index] = newData[i]
		finalNum++
	}
	return answer
}


func left(i int) int {
    return 2 * i+1
}

func right(i int) int {
    return 2*i + 2
}

func maxHeapify(a []int, i int) []int {
    l := left(i)
    r := right(i)
    var max int
    if l < len(a) && l >= 0 && a[l] < a[i] {
        max = l
    } else {
        max = i
    }
    if r < len(a) && r >= 0 && a[r] < a[max] {
        max = r
    }
    if max != i {
        a[i], a[max] = a[max], a[i]
        a = maxHeapify(a, max)
    }
    return a
}

func buildMaxHeap(a []int) []int {
    for i := len(a)/2 - 1; i >= 0; i-- {
        a = maxHeapify(a, i)
    }
    return a
}

func heapSort(a []int) []int {
    a = buildMaxHeap(a)
    size := len(a)
    for i := size - 1; i >= 1; i-- {
        a[0], a[i] = a[i], a[0]
        size--
        maxHeapify(a[:size], 0)
    }
    return a
}*/

/*
func Push(h Interface, x interface{}) {
    h.Push(x)        // call to Push defined on your custom type
    up(h, h.Len()-1) // **heapification**
}

// Pop removes the minimum element (according to Less) from the heap
// and returns it. The complexity is O(log(n)) where n = h.Len().
// It is equivalent to Remove(h, 0).
//
func Pop(h Interface) interface{} {
    n := h.Len() - 1
    h.Swap(0, n)
    down(h, 0, n) // **heapification**
    return h.Pop()
}*/
