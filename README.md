# NMP-Hydro: A C#-language based Noah-MP model
## Introduction
NMP-Hydro is a hydrological model based on Noah-MP but coded using the CSharp (C#) language. This model was crafted by creating a framework and accurately translating the original Noah-MP LSM code from WRF-Hydro 3.0, coupled with a Muskingum method-based river routing model(Liu et al., 2023). CSharp, recognized for its modern and object-oriented approach, is widely used for software development across various platforms, particularly on the Windows operating system. 
NMP-Hydro offers several advantages over the original Noah-MP and WRF-Hydro. Unlike the original version that requires compiling for each computer and predominantly relies on Unix-like systems, NMP-Hydro can seamlessly run on Windows systems supporting the Microsoft Dotnet Framework. The executable files, once compiled, can be easily packaged and distributed to other Windows computers, providing convenience for users less familiar with Unix-like operations. The utilization of the CSharp language facilitates advanced software programs for code visualization and analysis, enhancing user convenience for code reading and modification. The model's design aligns with the input datasets and settings in the 'namelist' file, ensuring compatibility with WRF-Hydro 3.0. Both the translated Noah-MP LSM simulation and the river routing simulation in NMP-Hydro support parallel execution on common personal computers.

## The Development process
### Translation of Noah-MP that embeded in WRF-Hydro 3.0
Our primary focus in developing NMP-Hydro involved translating the original FORTRAN code of Noah-MP into the C# language. The overarching objective of this translation was to create a hydrological model based on Noah-MP capable of functioning seamlessly on Windows systems. It is essential to note that this translation is based on a relatively older version of Noah-MP utilized in WRF-Hydro 3.0, as the process commenced before the release of Noah-MP 5.0 (He et al., 2023).
Converting FORTRAN code into CSharp is not straightforward due to significant differences in syntax between the two languages. The reconstruction of the model in the CSharp language follows a straightforward object-oriented design. While FORTRAN is traditionally a function-based language, the core Noah-MP module's functions, subroutines, and state variables are encapsulated as members within a class named GridCell (Fig. 1(a)). This class represents all Noah-MP behaviors within a grid box. The variable names, function definitions, data structures, and execution logic have been kept largely consistent with the original FORTRAN code, ensuring user-friendliness for those familiar with Noah-MP. To handle multiple grid boxes, another class named Driver is employed. This class manages tasks such as initializing model variables, creating multiple grid boxes, reading/writing files, and controlling the execution of the model.
Throughout the translation process, a key focus was addressing operations on FORTRAN arrays (Fig. 1(b)), crucial for representing the state of soil and snow layers in Noah-MP. Unlike CSharp, FORTRAN allows arrays to have user-specified index ranges (e.g., index values from -3 to 4). However, in CSharp, the first index of all arrays invariably starts from 0. To streamline the translation, we introduced a new array class named FortArray, designed to mimic FORTRAN arrays. The inner array data in FortArray adheres to standard CSharp conventions, accepting 0 as the inner index of the first element. Yet, externally, the class restricts access to the array values through additional indices. The class provides methods for index translation from extra indices to inner indices.
For instance, if a FORTRAN array of 8 elements has an index range from -3 to 4, this array is translated into a FortArray. The FortArray has a standard inner array of 8 elements, accompanied by two arguments representing the start index and the end index. For FortArray objects, encountering an extra index of -3 results in the FortArray class determining an inner index of 0 by subtracting the start index. Simultaneously, an extra index of 4 is translated to an inner index of 7. This array translation technique ensures that all the original execution logic in Noah-MP is seamlessly preserved in NMP-Hydro.

### Parallelization on PCs
The model boasts support for parallel execution, implemented through the native parallel functions of the CSharp language. These functions efficiently allocate computational tasks for distinct grid boxes to different CPU threads. For instance, if a specific domain requires the execution of 2400 grid boxes, and the tasks are assigned to 8 threads, each thread is responsible for completing the tasks of approximately 300 grid boxes. It's crucial to note that if the number of specified threads exceeds the actual number of CPU cores, multiple threads may end up executing on a single CPU core. Therefore, specifying more threads than the available CPU cores does not contribute to an overall improvement in execution speed. 

### Coupling with a river routing module
A parallel river routing module based on the Muskingum method (Liu et al., 2023) was added into the model, deviating from the previous utilization of the coupled RAPID model in WRF-Hydro (Lin et al., 2018). This parallel river routing module, implemented using CSharp, incorporates our unique techniques: 
(1) An array-based sequential processing method for Muskingum routing.
(2) A straightforward equal-sized domain decomposition method.
(3) Three distinct parallelization schemes for river routing.
(4) A specific sorting approach for river segments used in domain decomposition.
This approach's primary advantage lies in its ability to straightforwardly decompose any river network into multiple domains with an equal number of river segments. Achieved by evenly dividing the river segment list into any number of blocks, this innovation capitalizes on the inherent tree-like structure present in most river networks. Importantly, it does not necessitate consideration of the topological conditions specific to a given network, as required in studies such as Mizukami et al. (2021) or David et al. (2015). This design allows parallel execution of river routing on modern personal computers equipped with multi-CPU cores.
The integration of the river routing module with the Noah-MP LSM involves assigning lateral inflows from the LSM-simulated total runoff to the river routing model. In the present NMP-Hydro configuration, we utilize a straightforward catchment centroid-based coupling interface (David et al., 2015). This method designates the LSM grid cell containing the catchment centroid (referred to as the "centroid cell") as the location for a river reach to receive lateral inflows. At a specific temporal step, the computed contributing runoff discharge Qlat (unit: m3/s) is determined by the following expression: 

Q_lat=R(nx,ny)×F×1000


where R(nx,ny) is the runoff (mm, surface + subsurface) simulated by the LSM during the time step, F is the catchment area (km2) contributing water to the current river segment.
Alternatively, employing weighted assignments from different grid boxes, akin to the method utilized in (Lin et al., 2018), is also a valid approach. However, this method requires the generation of weights from multiple grid boxes. Given the size of each grid box (equivalent to the resolution of meteorological data, typically ranging from 25 km to 100 km), and considering that each grid box can encompass the catchment areas of multiple river segments, the coupling approach using area weighting is unlikely to yield substantial improvements for most river segments. 

## Usage
Csharp version Noah-MP configuration and usage instructions: 
We need to configure the wrfinput file for the simulation area, such as the wrfinput_YellowRiver_0.05_SGBP.nc file, which can be set in the namelist.Hrldas file (see figure below). This is a data file in netcdf format. 
To prepare the wrfinput file, it is necessary to use WPS's geogrid.exe to generate the geo_ em_ *. nc file, and then use a script with the code nco_stcript_to-create_rfinput_from_geogrid-with_nco_initification. sh to convert and generate it. 
The end of the file name Wrfpinput.d02 should indicate the IGBP or USGS vegetation classification code.

the setting of wrfinput file and the driving data folder:
![图片](https://github.com/user-attachments/assets/25e79c8c-8414-45cf-8628-8f65cf807a14)

Driver data requirements:

The driving data is meteorological data that includes 7 variables (wind temperature, humidity, pressure, precipitation, radiation). File every 3 hours. The file name contains time information, such as: 2003010112.LDASIN-DOMAIN1
Place the weather driver files (every three hours) in the directory pointed to by INDIR. These files can be placed directly in this directory or in a directory named after each year (such as directory 2000, 2001, 2002).

 ![图片](https://github.com/user-attachments/assets/bbab0a52-1b78-4d59-86d3-195210fdace7)

 
