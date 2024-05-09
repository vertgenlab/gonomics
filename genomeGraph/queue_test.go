package genomeGraph

import (
	"container/heap"
	"testing"
)

// Sample SeedDev for testing
var testSeed1 = SeedDev{TotalLength: 10}
var testSeed2 = SeedDev{TotalLength: 5}
var testSeed3 = SeedDev{TotalLength: 15}

// ------------------
// Tests for MaxItem 
// ------------------

// No specific tests for MaxItem itself, as it's primarily a struct for data holding

// ----------------------------
// Tests for MaxPriorityQueue 
// ----------------------------

func TestPushAndPop(t *testing.T) {
    pq := make(MaxPriorityQueue, 0)

    heap.Push(&pq, &MaxItem{Value: &testSeed3, Priority: 15})
    heap.Push(&pq, &MaxItem{Value: &testSeed1, Priority: 10})
    heap.Push(&pq, &MaxItem{Value: &testSeed2, Priority: 5})

    // Verify order
    item := heap.Pop(&pq).(*MaxItem)
    if item.Value != &testSeed3 || item.Priority != 15 {
        t.Errorf("Error: Incorrect top element. Expected %v, got %v", &testSeed3, item.Value)
    }

    item = heap.Pop(&pq).(*MaxItem)
    if item.Value != &testSeed1 || item.Priority != 10 {
        t.Errorf("Incorrect top element. Expected %v, got %v", &testSeed1, item.Value)
    }
}

func TestLen(t *testing.T) {
    pq := make(MaxPriorityQueue, 0)
    heap.Push(&pq, &MaxItem{Value: &testSeed1, Priority: 10})
    heap.Push(&pq, &MaxItem{Value: &testSeed2, Priority: 5})

    if pq.Len() != 2 {
        t.Errorf("Error: Incorrect length. Expected 2, got %d", pq.Len())
    }
}
