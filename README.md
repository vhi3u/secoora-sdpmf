# secoora-sdpmf
SECOORA Slocum Data Processing and Management Framework

This GitHub repository contains useful tools to support the [SECOORA Glider Observatory](https://secoora.org/data/secoora-glider-observatory/) in Slocum glider processing and management for consistent and reliable practice. 

SDPMF supports data processing and management for the following gliders across different generations:
* G1: pelagia, salacia, bass, sam, ramses (lost)
* G2: modena (lost)
* G3: angus, franklin, unit_1091


Processing and output files are organized based on the level of processing done. 

L0: "Bare minimum" work done to convert ASCII .ebdasc/.dbdasc files to .mat files for further processing. Glider data is typically stored in .dbd (flight) and .ebd (science) files that are dumped from the glider itself, which is then converted to ASCII. The process to pull .dbd and .ebd to ASCII (usually in the form of .dbdasc and .ebdasc) are not covered in this repo. 

L1: Essential corrections and processing from the L0.mat files. Processing is split across three sensors typically found in the Slocum gliders
* CTD (Conductivity, Temperature, Depth): Thermal lag correction based on [Garau et al., 2011](https://journals.ametsoc.org/view/journals/atot/28/9/jtech-d-10-05030_1.xml). 
* ECO (Chlorophyll, CDOM, Backscatter, etc.): Minmax filter
* DO (Dissolved Oxygen): [WIP]

L2: Gridding L1 data. Takes in all L1 (CTD, ECO, and DO) into one L2.mat. 

