
[![Go Report Card](https://goreportcard.com/badge/github.com/vertgenlab/gonomics)](https://goreportcard.com/report/github.com/vertgenlab/gonomics)
[![codecov](https://codecov.io/gh/vertgenlab/gonomics/branch/main/graph/badge.svg?token=SLasptsu7B)](https://codecov.io/gh/vertgenlab/gonomics)

# gonomics
A collection of genomics software tools written in Go (golang).

![gonomicsFigure](https://github.com/vertgenlab/gonomics/assets/49315918/7475cdf8-c20c-45e3-a0c1-3c2fbf6ef7a7)

### Complete Gonomics Documentation
The complete documentation for Gonomics can be found [here](https://pkg.go.dev/github.com/vertgenlab/gonomics).

### Installation

**1. [Install go.](https://go.dev/doc/install)**

**2. Install gonomics**

*Option 1: Executables only* - `go install github.com/vertgenlab/gonomics/...@latest`  

*Option 2: Complete library & executables*  
```
## Clone gonomics repository  
git clone https://github.com/vertgenlab/gonomics.git && cd gonomics

## Run gonomics tests
go test ./...

## Install executables
go install  ./...
```

Executables will be present in the Go binary folder (`~/go/bin` by default)
- Running a command with no arguments will print the usage statement for that command: `~/go/bin/commandName`

The command line tools' code is located in `/gonomics/cmd/`

* More information about using our tools can be found here: [Using gonomics commands](https://github.com/vertgenlab/gonomics/wiki/Using-gonomics-commands)
* To see every command in gonomics in command groups with their usage listed run `~/go/bin/gonomics`
* Instructions and examples on making your own tools can be found here: [Writing new tools with gonomics](https://github.com/vertgenlab/gonomics/wiki/Writing-new-tools-with-gonomics)
---

### To create a docker container of Gonomics

1. Inside the root directory of the gonomics repository run this command to build a docker container of gonomics:

```
docker build . -t gonomics
```

* This command takes a working directory as an input. A period is used here to indicate the current directory. Any scripts, data, or files you reference will build local paths from that directory.
* `-t` (tag) is used to name you container. As you acculate more docker images, you can easily reference them by their assgined tags.

2. Inside the container, you are essently running the same set up for build from source. I set up my working directory as: WORKDIR: `/src/github.com/vertgenlab/gonomics`

* alternatively you could just set it as `$HOME/src/github.com/vertgenlab/gonomics`, but this makes more sense to me because the general convention is to create your applications at the `/app` or `/gonomics` but as a go programing building in src is how I think about it.

3. Once we set our working directory in the container, we need to add code we wrote with `ADD . .` which takes a local directory as the first argument and the destination (inside the docker container) as the second. In other words we are recursively copying local files in gonomics into our `WORKDIR` which we set in step 2.

* `ADD (everything in gonomics dir) to here /src/github.com/vertgenlab/` is what that means.

4. I think everything should be familiar at this point you run the following commands to install all package dependencies with:

```
    go mod download
    go mod vendor
    go mod verify
```

5. Finally, `go install ./...` will build all gonomics executables into your $GOBIN which is located at `/go/bin` by default.

6. Now you can basically run any gonomics cmd globally in the container:

```
    docker run gonomics chimpAncestorRecon -h
```
    * or use the containers golang compiler globally on your desktop without having to worry about setting your $GOROOT, $GOPATH, $GOBIN

* all your gonomics dependencies should be all set up at this point

```

docker run -v $(pwd):/mnt gonomics go run $script.go

```

* this command will mount your local file system into the `/mnt` directory of your container.

---

### Compatibility with Previous Golang Versions
<p>Gonomics is compatible with Golang version 1.18 and above. Please note that due to changes in random number generation since Golang v1.20,
many cmd tests will fail if run on 1.18 or 1.19, including cmds that use random numbers, such as the simulate commands and MCMC sampling.
However, we expect these programs to function as intended in these older versions.</p>

---

### Authors:

* Eric Au
* Luke C. Bartelt
* Olivier Boivin
* Sophie Campione
* Christiana Fauci
* Craig B. Lowe
* Yanting Luo
* Riley J. Mangan
* Chelsea R. Shoben
* Daniel A. Snellings
* Seth Weaver
