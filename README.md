# QualMS
A prototype Mass Spectrometry Quality Analysis Tool for replicate MS data.
Designed by Scott Walmsley, Damian Fermin, and Hyungwon Choi.
Algorithmic design by Hyungwon Choi (BIC)
Informatic design and code by Damian Fermin.
Overall Design and code by Scott Walmsley.

#Usage
```
USAGE: qualMS.exe -k -b -d
```
```
	Command flags:<br>
    -k           Assume 'k' states (Run -b option first)<br>
    -d <path>    Path to folder containing input matrix files. ** Required input **<br>
    -b <#:#>     Provide an estimate for best number of 'k' states using BIC\<br>
    ```
    
    Example usage: -b 3:10 means consider between 3 and 10 states<br>
```
	     Using the -b flag, qualMS will only report the best value for K and exit<br>

```
