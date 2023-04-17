// Package exception provides globally available functions to simplify exit after unrecoverable errors are encountered.
package exception

import (
	"log"
	"math/rand"
	"os"
	"strings"
)

// init is run as the exception package is imported and modifies the GODEBUG variable to contain
// randautoseed=0 so that functions in gonomics that call the rand package are deterministic by default.
// The seed for the rand package can still be changed in individual packages or cmds by calling rand.Seed(num).
func init() {
	dbgVar := os.Getenv("GODEBUG")

	if dbgVar == "" {
		err := os.Setenv("GODEBUG", "randautoseed=0")
		PanicOnErr(err)
		return
	}

	variables := strings.Split(dbgVar, ",")
	var found bool
	for i := range variables {
		if strings.HasPrefix(variables[i], "randautoseed=") {
			variables[i] = "randautoseed=0"
			found = true
		}
	}

	if !found {
		variables = append(variables, "randautoseed=0")
	}

	err := os.Setenv("GODEBUG", strings.Join(variables, ","))
	PanicOnErr(err)
	rand.Seed(0)
}

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

// WarningOnErr will output a warning message (WARNING: err text) if err != nil.
func WarningOnErr(err error) {
	if err != nil {
		log.Printf("WARNING: %s", err)
	}
}
