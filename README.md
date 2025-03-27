```
╭━━━┳━━━┳━━━╮
┃╭━╮┃╭━╮┣╮╭╮┃
┃╰━━┫╰━━╮┃┃┃┣━━┳╮ ╭╮
╰━━╮┣━━╮┃┃┃┃┃╭╮┃┃ ┃┃
┃╰━╯┃╰━╯┣╯╰╯┃╰╯┃╰━╯┃
╰━━━┻━━━┻━━━┫╭━┻━╮╭╯
            ┃┃ ╭━╯┃
            ╰╯ ╰━━╯

SSDpy - Structural & Stochastic Dynamics
Tools for the structural dynamic and stochastic analysis of structures
V. Denoël. Toolbox for research and teaching.
```



## Stabilized versions
Date     ┃  Release  ┃  Comments
May  22  ┃  V.0.0.0  ┃  Development of the main folder arrangement
Mar. 23  ┃  V.0.0.1  ┃  Add dyn (integrators) and basic signal processing
Dec. 24  ┃  V.0.0.2  ┃  Add wind engineering tools, BAMM

##Contents
bayes  # bayesian updating
constitutive # constitutive laws
dyn    # structural dynamics : integrators, ...
gcdc   # interacting with GCDC sensors
misc   # miscelaneous
plot   # our plot formatting tools
sandbox # play here!
signal # some signal processing techniques : first passage time, ...
stat   # advanced statistics methods
tools  # miscellaneous tools for plotting, time-freq analysis
wind   # wind engineering : buffeting, turbulence simulator, flutter analysis, ...

## Installation

### For Regular Users
To install the package, follow these steps:
```bash
# Navigate to your target folder
cd your_target_folder

# Clone the repository
git clone https://github.com/vincentdenoel/SSD_tools.git

# Navigate into the cloned folder
cd SSD_tools

# Install the package
pip install .
```
### For developers and contributors:
```bash
# Clone the github project on a local drive
cd your_target_folder
git clone https://github.com/vincentdenoel/SSD_tools.git

# Move to your SSDpy environment:
conda activate SSDpy

# Install the package to your venv in editable mode:
pip install -e your_target_folder/SSDpy
```
