# hamFoam

An open-source solver for coupled heat and moisture transport in porous media based on OpenFOAM

hamFoam is mainly intended for the absorption, transport and storage simulations of coupled heat and moisture within porous building materials.

<img src="https://carmeliet.ethz.ch/research/downloads/coupled-heat-and-moisture-transport-solver-for-openfoam/_jcr_content/par/fullwidthimage/image.imageformat.fullwidth.1496181206.png"  width="400">

The solver is tested for the following OpenFOAM versions:

* OpenFOAM-org (OpenFOAM Foundation) v6, v7, v8, v9, v10, v11
* OpenFOAM-com (OpenCFD-ESI) v1806

### Usage

You can compile the solver for a specific OpenFOAM version by checking out the commit with corresponding tag. For example, for OpenFOAM v9:

	git clone https://github.com/OpenFOAM-BuildingPhysics/hamFoam.git
	cd hamfoam
	git checkout tags/of-org_v9.0
	wmake

See the list of tags for different versions [here](https://github.com/OpenFOAM-BuildingPhysics/hamFoam/tags)

### Documentation

Read [hamFoam wiki](https://gitlab.ethz.ch/openfoam-cbp/solvers/hamfoam/-/wikis/home) for documentation.

### Tutorial cases

Tutorial cases modeling HAMSTAD Benchmark case 4 (response analysis) and case 5 (capillary active inside insulation) can be found [here](https://github.com/OpenFOAM-BuildingPhysics/hamFoam-tutorials).

### Validation

The solver results have been compared with the benchmark cases from HAMSTAD (Heat, Air and Moisture StanDardization). Please see [tutorials](https://github.com/OpenFOAM-BuildingPhysics/hamFoam-tutorials) to download and run cases 4 and 5 from HAMSTAD.

<i>Hagentoft, C-E., 2002. HAMSTAD – Final report: Methodology of HAM-​modeling,
Report R-​02:8. Gothenburg, Department of Building Physics, Chalmers University
of Technology.</i>

More information at the Chair of Building Physics: https://carmeliet.ethz.ch
