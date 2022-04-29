# mpMRI_TP_DP

**Author**: Caterina Brighi

Scripts for combining mpMRI data to obtain tumour probability maps to use for dose painting.

**Setup/Build/Install** 
To be able to use the above scipts you need to have the following python packages intalled:
matplotlib
SimpleITK
numpy
pandas
datetime
os
glob
gzip
shutil
xlsxwriter
scipy
dicom2nifti
skimage

You also need to have the ImageAnalysisFunctions.py, ImageResultsVisualizationFunctions.py and ImageStatisticsFunctions.py files in the same folders as the other .py files of this repo. 

**Usage** 
The main project file of this repository, containing a pipeline to use functional MRI data to develop tumor probability maps (according to the method in this publication [Verburg N, Koopman T, Yaqub MM, et al. Improved detection of diffuse glioma infiltration with imaging combinations: a diagnostic accuracy study. Neuro Oncol. 2020;22(3):412-422. doi:10.1093/neuonc/noz180](https://academic.oup.com/neuro-oncology/article/22/3/412/5573318)) and dose painting prescriptions is the mpMRI_QIN.py file. All the other files perform other tasks associated with this project.

**Directory Structure** 
NA

**Citation**
If you use code from this repository, please cite [![DOI](https://zenodo.org/badge/380894935.svg)](https://zenodo.org/badge/latestdoi/380894935)
