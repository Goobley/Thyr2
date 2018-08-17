# Thyr

More complete documentation to follow, though the public API is already documented in the `doc` folder.
Running Thyr:
- Install [Torch7](www.torch.ch)
- Make sure you have a C++14 compatible compiler
- Modify `compile.sh` to work with your environment
- Run `compile.sh`
- Clone [Plyght](https://github.com/Goobley/Plyght) into Thyr's folder, making sure you have the Python pre-requisites (numpy, matplotlib, pyqt5 (Mac only)).
- Inspect and/or modify the `ThyrDriver.lua` file (should be well commented)
- Create a directory with the same name as the `prefix` variable in `ThyrDriver.lua` in the working
  directory (e.g. `c7DataHR` in the default case).
- Run Plyght (i.e. `python3 Plyght/plyght.py`) and leave in running in background
- Run the simulation `th ThyrDriver.lua` (after whatever variant of `torch_activate` is required on your system)
- Post-processing of the simulations from the backed up "datacubes" is done using `th ThyrPostProcess.lua`.
  Set `prefix` in `ThyrPostProcess.lua` to the same string as was used in `ThyrDriver.lua`.

**NB: The default simulation in ThyrDriver.lua is currently very high resolution and will take a long time to compute. To get quick results from Thyr I would recommend `numVox = 32` with a refinement factor of 2 in `Thyr.create_high_res_footpoints`. To obtain results quickly it is also a good idea to reduce the number of frequencies being computed, especially the higher frequencies.

All files under [MIT License](https://opensource.org/licenses/MIT) (c) 2015-2018 Christopher Osborne, University of Glasgow, with the exception of `maf.lua` (c) Bjorn Swenson (MIT License) and the excellent [Penlight](http://stevedonovan.github.com/Penlight) libraries (c) Steve Donovan (MIT License).