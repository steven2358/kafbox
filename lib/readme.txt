How to format an algorithm
--------------------------

1. Make a copy of lib/dummy.m and name it MY_ALGORITHM.m
2. Open the file and replace all mentions of "dummy" to MY_ALGORITHM
3. In the first properties section, fill in the parameter values used by the algorithm and some default values
4. In the second properties section, fill in the variables that will be calculated by the algorithm
5. [optional] Adjust the "evaluate" method
6. Add the code for the algorithm to the "train" method

Check if the code is formatted correctly by running the script unit_test('MY_ALGORITHM') in the folder lib/test/

Some information on the code structure:
- Each algorithm is contained in a single file in the /lib folder. 
- Each algorithm is implemented as objects in matlab, using the classdef syntax.
- The object code contains only one iteration of the algorithm. The for-loop over the time index that governs the online operation goes in an external script.


How to contribute an algorithm to the toolbox
---------------------------------------------

Option 1: email it to me (steven@gtas.dicom.unican.es)

Option 2: fork the toolbox on GitHub (https://github.com/steven2358/kafbox), push your change to a named branch, then send me a pull request.