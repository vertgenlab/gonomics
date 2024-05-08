package interval

import (
	"encoding/gob"
	"io"
)

// WriteTree gobs the input tree and writes to w
func WriteTree(w io.Writer, tree Tree) error {
	comprehensiveRegister(tree)
	encoder := gob.NewEncoder(w)
	err := encoder.Encode(tree)
	if err != nil {
		return err
	}
	return nil
}

// ReadTree reads a gobbed tree from r to the returned Tree struct
func ReadTree(r io.Reader) (Tree, error) {
	decoder := gob.NewDecoder(r)
	var ans Tree
	err := decoder.Decode(&ans)
	if err != nil {
		return nil, err
	}
	return ans, nil
}

// comprehensiveRegister registers the type of EVERY data record in the tree.
// This is only necessary when the try contains multiple data types.
// NOTE: benchmarking shows that registration is not very computationally expensive
// so we will stick to comprehensiveRegister for the time being.
func comprehensiveRegister(tree Tree) {
	for _, val := range tree {
		for i := range val.Data {
			gob.Register(val.Data[i])
		}
	}
}

// easyRegister only registers the type of the first data record in each chromosomes tree.
// This is effective when the tree is composed of only a single data type.
func easyRegister(tree Tree) {
	for _, val := range tree {
		gob.Register(val.Data[0])
	}
}
