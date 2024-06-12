package genomeGraph

type SeedHeap struct {
	*Seed
}

type PriorityQueue []*SeedHeap

func (pq PriorityQueue) Len() int           { return len(pq) }
func (pq PriorityQueue) Less(i, j int) bool { return pq[i].TotalLength > pq[j].TotalLength }
func (pq PriorityQueue) Swap(i, j int)      { pq[i], pq[j] = pq[j], pq[i] }

func (pq *PriorityQueue) Push(x interface{}) {
	item := x.(*SeedHeap)
	*pq = append(*pq, item)
}

func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	*pq = old[0 : n-1]
	return item
}

func (pq PriorityQueue) Top() interface{} {
	if len(pq) == 0 {
		return nil
	}
	return pq[0] // The first element is the highest priority in a max-heap
}
