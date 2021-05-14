// Package exception provides globally available functions to simplify exit after unrecoverable errors are encountered.
package exception

import "log"

// PanicOnErr will panic if input error is not nil.
// Panic is the preferred kill method for packages in the codebase (i.e. everything not in the cmd folder).
func PanicOnErr(err error) {
	if err != nil {
		log.Panic(err)
	}
}

// FatalOnErr will fatal if input error is not nil.
// Fatal should generally be reserved for end user functions (i.e. in the cmd folder).
func FatalOnErr(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

// WarningOnErr will output a warning message (WARNING: err text) if err != nil
func WarningOnErr(err error) {
	if err != nil {
		log.Printf("WARNING: %s", err)
	}
}
