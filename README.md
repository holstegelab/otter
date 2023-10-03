# *otter*

## Dependencies

* A c++17 supporting compiler
* [htslib-1.3+](https://www.htslib.org/download/)
* [WFA2-lib](https://github.com/smarco/WFA2-lib) (included during install)
* [hclust-cpp](https://github.com/cdalitz/hclust-cpp/tree/master) (included)

Note that *htslib* needs be available in your library paths. You can do this by adding the *htslib* path in your *.bashrc* or *.bash_profile* file:

```bash
export CPATH=$CPATH:<your_local_htslib_path>
export LIBRARY_PATH=$LIBRARY_PATH:<your_local_htslib_path>
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<your_local_htslib_path>
```

## Installation

```bash
git clone https://github.com/holstegelab/otter.git && cd otter
mkdir build
make packages
make
```

## Usage

See commands and version by running:
```bash
./build/otter
```