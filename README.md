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