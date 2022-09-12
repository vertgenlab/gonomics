
[![Go Report Card](https://goreportcard.com/badge/github.com/vertgenlab/gonomics)](https://goreportcard.com/report/github.com/vertgenlab/gonomics)
[![codecov](https://codecov.io/gh/vertgenlab/gonomics/branch/main/graph/badge.svg?token=SLasptsu7B)](https://codecov.io/gh/vertgenlab/gonomics)

# gonomics

A collection of genomics software tools written in Go (golang).
Authors:

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

### Instructions for installing Golang (Go)

<p>Instructions are written for bash (common for Linux/masOS). Windows users can follow directions on this website. https://golangdocs.com/install-go-windows </p>
<br>

**1. Download Go**

<ul>
    Determine where you want the go program to install. <br>
    This will be separate from your directory for storing and working on go code (aka your go workspace). <br>
    Recommended: /usr/local <br>

    cd /pathWhereYouWantGoInstalled/
    wget <link for go version of choice> 

</ul>

[Versions of Go can be downloaded from here](https://golang.org/dl/)   - we are currently using 1.18

<ul>

    tar -xzf [version of go you downloaded]
</ul>
<br>
<br>

**2. Clone gonomics into Go/set up Go workspace**

*Gonomics installation with version control (for both users and contributors of gonomics):*
<ul>

    
    cd /pathWhereYouWantGonomicsCloned/
    mkdir -p src/github.com/vertgenlab
    cd src/github.com/vertgenlab
    git clone https://github.com/vertgenlab/gonomics.git
<ul>
This will download the repository to your current directory
</ul>

    cd gonomics
    go test ./…
<ul>
This should print to screen a line for each test with "ok" printed in the left margin when something passes.
</ul>

    go install ./...
<ul>
This will tidy up the necessary modules
</ul>
</ul>
<br>
<br>

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
