# Feature-preserving Isosurface Extraction


<p align="center"> 
<img src="teaser.png" width="1000">
</p>

This code is accompanied with the following work for extracting the iso surface of a given implicit surface:

Hao-Xiang Guo, Yang Liu, Hao Pan, Baining Guo. <b>Implicit Conversion of Manifold B-Rep Solids by Neural Halfspace Representation</b>, <i>ACM Transactions on Graphics(SIGGRAPH ASIA)</i>, 2022. 

[[Paper]](https://arxiv.org/abs/2209.10191)  [[Project Page]](https://guohaoxiang.github.io/projects/nhrep.html)


## Installation Guide （64-bits system) ##
The code was tested under both Linux (Ubuntu 18.04, 20.04) and Windows (Windows 10+), with and without GPU support (Cuda 10+). The following guide is for Linux system with Cuda GPU support. 

* Step 1: download essential packages [LibTorch 1.7.1+cu101](https://download.pytorch.org/libtorch/cu110/libtorch-cxx11-abi-shared-with-deps-1.7.1%2Bcu110.zip) and [CuDNN 8.0.5+cu101](https://developer.nvidia.com/rdp/cudnn-archive）. Extract the above packages into two folders `libtorch` and `cudnn` under the folder `IsoSurfacing`. 

    <i>Higher versions of LibTorch, CuDNN and CUDA can be used here, but you may need to install additional packages like `intel-mkl`, and specify CMake path.</i>     

* Step 2: Compile the code. You need `cmake` (ver >= 3.9) and `g++` (>= 8.0).
```bash
    cd IsoSurfacing && mkdir build && ./build.sh
```

## Usage:  Isosurface extraction from neural implicits

Run the following command to see the guidance.
```bash
    cd .IsoSurfacing/build/App/console-pytorch
    ./ISG_console_pytorch -h
```
The following message is printed on screen.
```
    --- Usage:
    ISG_pytorch [OPTION...]

    -i, --input arg      input model(*.pt)
    -d, --depth arg      octree depth (2-10)(default: 7)
    -b, --box arg        bounding box size (default: 2)
    -a, --angle arg      feature angle threshold for EMC  (default: 30 degree)
    -m, --method arg     Method: EMC, DC, MC (default: EMC)
    -o, --output arg     output mesh (ply format)
    -t, --threshold arg  threshold value for computing intersection. (default:
                        1e-7)
    -n, --maxiter arg    max iteration number for computing intersection.
                        (default: 50)
    -g, --gpu            use GPU model (default: true)
    -s, --setbuf arg     set buffer size (default: 131072)
    -x, --type arg       set model type (2dim or 3dim) (default: 2)
    -l, --alldepth       generate meshes for each depth layer (default, false)
    -v, --isovalue arg   isovalue (default: 0)
    -k, --vtkoput        output vtk format
    -c, --useoctree      use octree speedup (default, false)
    -y, --verbose        print progress (default, true)
    -h, --help           Print help
        --compare arg    GT mesh for SDF comparison (*.ply)
        --compres arg    resolution for for SDF comparison (default: 32)
```

Assume that you have a network model in pytorch (*.pt) format that takes a set of 3D coordinates as input and outputs the corresponding implicit values or signed distances,  you can use ``ISG_console_pytorch`` to extract the isosurface of the implicit function, with sharp feature preserving property (see our paper).

```bash
./ISG_console_pytorch -i ./test/teaser.pt -o outputmesh.ply -v -0.01 -d 8
```
Here ``-0.01`` is the specified isovalue, `-d 8` is the max depth of grid volume, you may need to set it to a high value if the isosurface contains very narrow regions that needs high-resolution grids. You can set `--useoctree=true` to speed up the computation, but it may fail to recover small regions in some circumstances occasionally. If you want to recover very flat feature curves, you can use `-a` and specify a small angle like `5` degree, if `EMC` method is employed.  You can set a large buffer size via `-s integervalue` if you own a GPU with much large memory. The default value `131072` is tested on a GPU with 11GB Mem.

Here `EMC` is the algorithm proposed by Kobbelt et al. --- [Feature Sensitive Surface Extraction from Volume Data, SIGGRAPH 2001](https://www.graphics.rwth-aachen.de/media/papers/feature1.pdf), `DC` corresponds to the Dual Contour algorithm proposed by Ju et al. --- [Dual Contouring of Hermite Data, SIGGRAPH 2022](https://people.engr.tamu.edu/schaefer/research/dualcontour.pdf).  The authors' codes are adopted and modified.

The program will output an isosurface mesh -- `outputmesh.ply`. It also generates `outputmesh_sharpedge.obj` (sharp edges detected by checking dihedral angles) and `outputmesh.mlp` (MeshLab project file). You can open `outputmesh.mlp` from MeshLab to visualize both the mesh and sharp edges.

## Usage:  Isosurface extraction for customized implicit functions

Please check the small project `App\console` and modify `App\console\MyImplicitFunc.h` to define your implicit function.  You only need to implement `is_inside`, `scalar_value` and `gradient` functions.
I also integrated other marching cube variants for testing, including `MC` -- classic marching cubes, `MC_33` -- Marching Cubes 33, `DMC` -- Dual Marching Cubes.
Note that `App\console` project uses CPUs only.

## Contact

Please contact me (Yang Liu yangliu@microsoft.com) if you have any question about this implementation. 

## Third Party Code in use

- [*IsoEx*](https://www.graphics.rwth-aachen.de/IsoEx/): src code of the paper --- [Feature Sensitive Surface Extraction from Volume Data, SIGGRAPH 2001](https://www.graphics.rwth-aachen.de/media/papers/feature1.pdf).
- [*OpenMesh*](https://www.graphics.rwth-aachen.de/software/openmesh/): OpenMesh Lib used by IsoEx.
- [*Dual Contouring*](https://sourceforge.net/projects/dualcontouring/): src code of the paper --- [Dual Contouring on Hermite Data, SIGGRAPH 2001](https://people.engr.tamu.edu/schaefer/research/dualcontour.pdf).
- [*cxxopts*](https://github.com/jarro2783/cxxopts): c++ command line option parser 
- [*termcolor*](https://github.com/ikalnytskyi/termcolor): print colored messages
- [*Happly*](https://github.com/nmwsharp/happly):  PLY file format parser.
- [*Dual MC*](https://github.com/dominikwodniok/dualmc):  The implementation of Dual Marching Cubes by Dominik Wodniok.
- [*TriangleMeshDistance*](https://github.com/InteractiveComputerGraphics/TriangleMeshDistance): signed distance function to a triangle mesh, used for evaluating isosurface quality.
- `Marching Cubes 33`: Marching Cubes 33 algorithm.
- `CIsoSurface`: An implementation of Marching Cubes algorithm.

