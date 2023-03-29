package main

import (
	"testing"
)

func TestGonomics(t *testing.T) {
	binPath, _ := getBin()
	buildCmdCache(binPath)
	getCache()
}
