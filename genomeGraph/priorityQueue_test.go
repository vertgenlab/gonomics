package genomeGraph

import (
	"container/heap"
	"testing"
)

// Sample Seed for testing
var testSeed1 = Seed{TotalLength: 10}
var testSeed2 = Seed{TotalLength: 5}
var testSeed3 = Seed{TotalLength: 15}

// ------------------
// Tests for SeedHeap
// ------------------

// No specific tests for SeedHeap itself, as it's primarily a struct for data holding

// ----------------------------
// Tests for PriorityQueue
// ----------------------------

func TestPushAndPop(t *testing.T) {
	pq := make(PriorityQueue, 0)

	heap.Push(&pq, &SeedHeap{Seed: &testSeed3})
	heap.Push(&pq, &SeedHeap{Seed: &testSeed1})
	heap.Push(&pq, &SeedHeap{Seed: &testSeed2})

	// Verify order
	item := heap.Pop(&pq).(*SeedHeap)
	if item.Seed != &testSeed3 || item.TotalLength != 15 {
		t.Errorf("Error: Incorrect top element. Expected %v, got %v", &testSeed3, item.Seed)
	}

	item = heap.Pop(&pq).(*SeedHeap)
	if item.Seed != &testSeed1 || item.TotalLength != 10 {
		t.Errorf("Incorrect top element. Expected %v, got %v", &testSeed1, item.Seed)
	}
}

func TestLen(t *testing.T) {
	pq := make(PriorityQueue, 0)
	heap.Push(&pq, &SeedHeap{Seed: &testSeed1})
	heap.Push(&pq, &SeedHeap{Seed: &testSeed2})

	if pq.Len() != 2 {
		t.Errorf("Error: Incorrect length. Expected 2, got %d", pq.Len())
	}
}
