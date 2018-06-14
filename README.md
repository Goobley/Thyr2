# Thyr

More complete documentation to follow, though the public API is already documented in the `doc` folder.
Running Thyr:
- Install [Torch7](www.torch.ch)
- Make sure you have a C++14 compatible compiler
- Modify `compile.sh` to work with your environment
- Run `compile.sh`
- Clone [Plyght](https://github.com/Goobley/Plyght) into Thyr's folder, making sure you have the Python pre-requisites (numpy, matplotlib, pyqt5 (Mac only)).
- Inspect and/or modify the `ThyrDriver.lua` file (should be well commented)
- Run Plyght (i.e. `python3 Plyght/plyght.py`) and leave in running in background
- Run the simulation `th ThyrDriver.lua` (after whatever variant of `torch_activate` is required on your system)


All files under [MIT License](https://opensource.org/licenses/MIT) (c) 2015-2017 Christopher Osborne, University of Glasgow, with the exception of `maf.lua` (c) Bjorn Swenson (MIT License) and `gauss_legendre.*` (c) 2007-2010 Pavel Holoborodko under modified BSD 3-clause license with additional no commercial use clause (see files for more information).
The excellent [Penlight](http://stevedonovan.github.com/Penlight) libraries are also used.