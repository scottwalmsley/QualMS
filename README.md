# QualMS
A prototype Mass Spectrometry Quality Analysis Tool for replicate MS data.
Designed by Scott Walmsley, Damian Fermin, and Hyungwon Choi.<br>
Algorithmic design by Hyungwon Choi (BIC)<br>
Informatic design and code by Damian Fermin.<br>
Overall Design and code by Scott Walmsley.<br>

#Usage
```
USAGE: qualMS.exe -k -b -d
```
```
	Command flags:<br>
    -k           Assume 'k' states (Run -b option first)
    -d <path>    Path to folder containing input matrix files. ** Required input **
    -b <#:#>     Provide an estimate for best number of 'k' states using BIC
    
```
```    
    Example usage: -b 3:10 means consider between 3 and 10 states
```
```
         Using the -b flag, qualMS will only report the best value for K and exit.

```
