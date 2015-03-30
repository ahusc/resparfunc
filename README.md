# resparfunc
Sage Code for the Restricted Partition Function

This is a sagemath implementation of an algorithm to construct the restricted partition function for a given list of positive numbers (also known als Sylvester's denumerant).

The code allows to explore the periodic structure of the partition function as described in the paper "Properties of the Restricted Partition Function" and implements the algorithm outlined therein. 

It comes in two versions - a simple one which is intended mainly for illustration of the algorithm outlined in the paper and an optimized version that can be used for restriction lists with more and / or greater numbers.

# Authors

The code has been written by Anne Huschitt based on ideas of Alfred Seiler. These ideas are described in the paper "Properties of the Restricted Partition Function" (to be published). 

# License

The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation. 

# Installation

A prerequisite for using this code is a properly working installation of [sage] (http://www.sagemath.org). After downloading the code, you may import it in a sage session by typing:

```python
sage: load("resparfunc.sage")
```
or, for the simple version:
```python
sage: load("simple-resparfunc.sage")
```

The optimized and the simple version cannot be used together in the same sage session.

# Usage

After importing the code, type
```python
sage: RestrictedPartitionFunction?
```
to get information about the usage of the package, including various examples. The docmentation can also be found on the [project website] (http://ahusc.github.io/resparfunc).

The simple version may not be usable for more and / or greater numbers. Use the
optimized version instead.

