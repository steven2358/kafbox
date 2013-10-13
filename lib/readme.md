Algorithm code structure
------------------------
- Each algorithm is contained in a single file in the `/lib` folder. 
- Each algorithm is implemented as an object in matlab, using the `classdef` syntax.
- Each algorithm has two basic public methods: `evaluate` and `train`.
- The object code contains only one iteration of the algorithm. The for-loop over the time index that governs the online operation should go in an external script.


Quick steps to code an algorithm in the toolbox' format
-------------------------------------------------------
1. Make a copy of `lib/kafbox_template.m` and name it `MY_ALGORITHM.m`.
2. Open the file and replace all mentions of "kafbox_template" to MY_ALGORITHM.
3. Clean up the file header by replacing all text between square brackets.
4. In the first properties section, fill in the parameter values used by the algorithm with their corresponding default values. Default values should be chosen such that the algorithm performs well on a whitened input signal (see `demo/demo_sinc_all.m`).
5. In the second properties section, fill in the variables that will be calculated by the algorithm.
6. [optional] Adjust the `evaluate` method.
7. Add the code for training the algorithm to the `train` method. It may be useful to rely on "helper" functions to keep the code modular. For coding style, see below.
8. Check if the code is formatted correctly by running the script `unit_test('MY_ALGORITHM')`.


How to contribute an algorithm to the toolbox
---------------------------------------------
Option 1: email it to me (steven@gtas.dicom.unican.es).

Option 2: fork the toolbox on GitHub (https://github.com/steven2358/kafbox), push your change to a named branch, then send me a pull request.


Coding style
------------
**Code** should be  
1. as human-readable as possible, in the first place;  
2. short and structured, in the second place.  

**Algorithm structure** should follow the pseudocode from the corresponding publication.

**Variable naming** should correspond to the nomenclature used in the corresponding publication whenever possible.

**Comments** should be used sparingly. Document the design and purpose of the code rather than its mechanics.
