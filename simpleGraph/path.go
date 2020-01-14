package simpleGraph

import (
	"fmt"
)

func AddPath(newPath uint32, allPaths []uint32) []uint32 {
	if allPaths == nil {
		allPaths = append(allPaths, newPath)
	} else if allPaths[len(allPaths)-1] == newPath {
		return allPaths
	} else {
		allPaths = append(allPaths, newPath)
	}

	return allPaths
}

func CatPaths(currPaths []uint32, newPaths []uint32) []uint32 {
	if len(newPaths) == 0 {
		return currPaths
	} else if len(currPaths) == 0 {
		return newPaths
	} else {
		currPaths = AddPath(newPaths[0], currPaths)
		currPaths = append(currPaths, newPaths[1:]...)
		return currPaths
	}
}

func reversePath(alpha []uint32) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func PathToString(allPaths []uint32, gg *SimpleGraph) string {
	var s string = ""
	//fmt.Printf("length of paths %d\n", len(allPaths))
	if allPaths == nil {
		return s
	} else {
		s += fmt.Sprint(gg.Nodes[allPaths[0]].Id)
		if len(allPaths) > 1 {
			for i := 1; i < len(allPaths); i++ {
				s += ":" + fmt.Sprint(gg.Nodes[allPaths[i]].Id)
			}
		}
	}
	return s
}
