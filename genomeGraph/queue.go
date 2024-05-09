package genomeGraph

type MaxItem struct {
	Value    *SeedDev
	Priority uint32
}

type MaxPriorityQueue []*MaxItem

func (pq MaxPriorityQueue) Len() int           { return len(pq) }
func (pq MaxPriorityQueue) Less(i, j int) bool { return pq[i].Priority > pq[j].Priority }
func (pq MaxPriorityQueue) Swap(i, j int)      { pq[i], pq[j] = pq[j], pq[i] }

func (pq *MaxPriorityQueue) Push(x interface{}) {
	item := x.(*MaxItem)
	*pq = append(*pq, item)
}

func (pq *MaxPriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	*pq = old[0 : n-1]
	return item
}
