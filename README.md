The benchmarks can be found in the example_* directories.
The main *c or *cpp files are incorporate with approximation in both compute and memory. 
Ideal approximation is obtained by gradient descent which is incorporated into these applications.
The gradient descent tuning is hardware cognizant - the applications call gem5 (a hardware simulator) which executes a snippet of the application in hardware and collect stats such as ipc and related metrics.
The error is modelled within the application itself via considering compute/memory errors.
Finally gradient descent is used to optimize the error-vs-ipc metric - so that the approximation is tuned for maximum increase in ipc for minimum increase in error.
Most directories contain RUN/COMPILE options to provide details for running and compiling these benchmarks.

Please refer / cite SHASTA (TACO '20) - https://dl.acm.org/doi/10.1145/3412375

 
