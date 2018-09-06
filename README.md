Copyright (C) 2018 King Abdullah University Of Science and Technology 

[High-Performance Visualization Group](http://vccvisualization.org/)
# Time-Dependent Flow seen through Approximate Observer Killing Fields

This directory contains the matlab code and test data to reproduce results presented in the work of Hadwiger et al. [1].

The code was tested with Matlab 2017b.

The script 'scriptOptimizeObserverField.m' computes an approximate observer Killing field u for a given input flow field v.

The script 'scriptResampleObservedField.m' resamples an input flow field v as it is seen through the approximate observer Killing field u, producing the observed field w_r for a given observation time r.

The output fields u and w_r are saved in raw binary format; float32 x,y,t ordered (x fastest).

# Running the Examples

To run the different experiments from the paper change the variable 'example' in scriptOptimizeObserverField and scriptResampleObservedField.
The examples 'Cylinder2D' and 'ocean' require additional datasets. [Download the data set archive](https://www.dropbox.com/s/t4k5te6usafailt/data.zip?dl=0 "data"),
extract the archive and put it in the 'data' folder before running the examples.

# Acknowledgements

The ocean data set is produced by SSALTO/DUACS, distributed by [AVISO, with support from CNES](http://www.aviso.oceanobs.com/duacs), and made available in the work of Haller et al. [2].
The Cylinder2D data set was simulated by Tino Weinkauf [3], using the Free Software 'Gerris Flow Solver' [4].

# Citation
```If you want to cite our work:
@article{Hadwiger2019ObserverKillingFields,
 title = {Time-Dependent Flow seen through Approximate Observer Killing Fields},
 author = {Hadwiger, Markus and Mlejnek, Matej and Theu{\ss}l, Thomas and Rautek, Peter},
 journal = {IEEE Transactions on Visualization and Computer Graphics (Proceedings IEEE Scientific Visualization 2018)},
 year = {2019}
 volume = {25},
 number = {1},
 pages = {to appear}
}
```

# Contact: 

Peter Rautek peter.rautek@kaust.edu.sa

Matej Mlejnek matej.mlejnek@kaust.edu.sa

# References:

[1] Hadwiger, M., Mlejnek, M.,  Theussl, T., Rautek, P., [Time-Dependent Flow seen through Approximate Observer Killing Fields](http://vccvisualization.org/research/killingobservers/). IEEE Transactions on Visualization and Computer Graphics (Proceedings IEEE Scientific Visualization 2018), 25(1), 2019.

[2] Haller, G., Hadjighasem, A., Farazmand, M., Huhn, F., Defining coherent vortices objectively from the vorticity. Journal of Fluid Mechanics, 795, 136-173, 2016. 

[3] Weinkauf, T., Theisel H., Streak Lines as Tangent Curves of a Derived Vector Field. IEEE Transactions on Visualization and Computer Graphics (Proceedings Visualization 2010), 16(6), 1225-1234, 2010.

[4] Popinet, S., Free Computational Fluid Dynamics. ClusterWorld, 2(6), 2004.