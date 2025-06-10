[![DOI](https://zenodo.org/badge/978270786.svg)](https://doi.org/10.5281/zenodo.15351035)

NucDraw is a simple Python package inspired by Forgi for the creation of NUCleic acid structures DRAWings.
It relies on viennaRNA and matplotlib to convert 1D dot-bracket structures into easily to interpret 2D drawings.

This package is well-suited to generate many structures in an automated way from multiplexed and high-throughput data, and allows to visualize multiple strands without pseudo-knots.

To install simply run from terminal:
- pip install nucdraw

If you use this package for your project, please remember to cite:
- 10.5281/zenodo.15352138

Several functionalities are included to allow the customization of your graphs.
See the following examples.


```python
from nucdraw import NucDraw
```


```python
# Let's generate a simple-to-read graph for a long RNA fold

seq = "CACAAUGUGGCCGAGGACUUUGAUUGCACAUUGUUGUUUUUUUAAUAGUCAUUCCAAAUAUGAGAUGCGUUGUUACAGGAAGUCCCUUGCCAUCCUAAAAGCCACCCCACUUCUCUCUAAGGAGAAUGGCCCAGUCCUCUCCCAAGUCCACACAGGGGAGGUGAUAGCAUUGCUUUCGUGUAAAUUAUGUAAUGCAAAAUUUUUUUAAUCUUCGCCUUAAUACUUUUUUAUUUUGUUUUAUUUUGAAUGAUGAGCCUUCGUGCCCCCCCUUCCCCCUUUUUUGUCCCCCAACUUGAGAUG"
mfe = ".((((((((....((....)).....)))))))).......................................((((..(((............................((((((....))))))...........((((((............))))))............)))..))))......................................................................................................................"

nc = NucDraw(mfe)
nc.generate(degree=90)
nc.plotter(8, bckwargs={'lw':2, 'color':'k'}, bpkwargs={'lw':2, 'c':'red'}, scwargs={'s':10, 'c':'k'})
```


    
![png](README_files/README_2_0.png)
    



```python
# Let's focus on one section and increase the details

seq = "CACAAUGUGGCCGAGGACUUUGAUUGCACAUUGUUGUUUUUUUAAUAGUCAUUCCAAAUAUGAGAUGCGUUGUUACAGGAAGUCCCUUGCCAUCCUAAAAGCCACCCCACUUCUCUCUAAGGAGAAUGGCCCAGUCCUCUCCCAAGUCCACACAGGGGAGGUGAUAGCAUUGCUUUCGUGUAAAUUAUGUAAUGCAAAAUUUUUUUAAUCUUCGCCUUAAUACUUUUUUAUUUUGUUUUAUUUUGAAUGAUGAGCCUUCGUGCCCCCCCUUCCCCCUUUUUUGUCCCCCAACUUGAGAUG"
mfe = ".((((((((....((....)).....)))))))).......................................((((..(((............................((((((....))))))...........((((((............))))))............)))..))))......................................................................................................................"

seq = seq[:50]
mfe = mfe[:50]

nc = NucDraw(mfe)
nc.generate(degree=90)
nc.plotter(8, bckwargs={'lw':2, 'color':'k'}, bpkwargs={'lw':2, 'c':'red'}, scwargs={'s':10, 'c':'k'})
nc.plot_circles(circle_size = 4, circle_color='white')
nc.plot_sequence(seq, {'fontsize':8, 'color':'k'})

# Let's annotate with numbers.
# nc.numbering_outside(distance, shape, number spacing, {'fontsize':10, 'color':'k'})
# "distance" determines how far the numbers will be from the structure
# "shape" tunes the direction along which the numbers will move away from the structure
# "number spacing" can either be an integer or a list.
# An integer "i" will cause every "i-th" number to be displayed.
# A list will determine exactly what numbers to be shown.

# We can add a number at a distance=15, with shape=5 and spacing=10

nc.numbering_outside(15, 5, 10, {'fontsize':10, 'color':'k'})

```


    
![png](README_files/README_3_0.png)
    



```python
# Let's focus on one section and increase the details
# Let's color-code the nucleobases

seq = "CACAAUGUGGCCGAGGACUUUGAUUGCACAUUGUUGUUUUUUUAAUAGUCAUUCCAAAUAUGAGAUGCGUUGUUACAGGAAGUCCCUUGCCAUCCUAAAAGCCACCCCACUUCUCUCUAAGGAGAAUGGCCCAGUCCUCUCCCAAGUCCACACAGGGGAGGUGAUAGCAUUGCUUUCGUGUAAAUUAUGUAAUGCAAAAUUUUUUUAAUCUUCGCCUUAAUACUUUUUUAUUUUGUUUUAUUUUGAAUGAUGAGCCUUCGUGCCCCCCCUUCCCCCUUUUUUGUCCCCCAACUUGAGAUG"
mfe = ".((((((((....((....)).....)))))))).......................................((((..(((............................((((((....))))))...........((((((............))))))............)))..))))......................................................................................................................"

seq = seq[:50]
mfe = mfe[:50]

nc = NucDraw(mfe)
nc.generate(degree=90)
nc.plotter(8, bckwargs={'lw':2, 'color':'k'}, bpkwargs={'lw':2, 'c':'k'}, scwargs={'s':10, 'c':'k'})
nc.plot_circles(seq, circle_size = 4)
nc.plot_sequence(seq, {'fontsize':8, 'color':'k'})
```


    
![png](README_files/README_4_0.png)
    



```python
# Let's draw a 2-strands complex

seq1 = 'UGACGUAAAACUGAC'
seq2 = 'UGUUACCGUA'
seq = "".join([seq1, seq2])
mfe = '..((((..(((....+.))))).)).'

nc = NucDraw(mfe)
nc.generate()
nc.plotter(6, bckwargs={'lw':2, 'color':'k'}, bpkwargs={'lw':3, 'c':'k'}, scwargs={'s':10, 'c':'k'})
nc.plot_circles(seq, circle_size = 3, circle_color='white')
```


    
![png](README_files/README_5_0.png)
    



```python
# Let's draw a 3-strands complex and color the strands differently

seq1 = 'UGACGUAAAACUGAC'
seq2 = 'UGUUACCGUAGUACG'
seq3 = 'ACCGUAC'
seq = "".join([seq1, seq2, seq3])
mfe = '..((((..(((....+.))))).)).(((((+..)))))'

nc = NucDraw(mfe)
nc.generate()
nc.plotter(8, bckwargs={'lw':2, 'color':'k'}, bpkwargs={'lw':3, 'c':'k'}, scwargs={'s':10, 'c':'k'})
nc.plot_circles(seq, circle_size = 2, circle_color='white')
nc.multistrand_coloring(clr=['red', 'blue', 'green'], bckwargs={'lw' : 3})

# Let's annotate with numbers.
# nc.numbering_outside(distance, shape, number spacing, {'fontsize':10, 'color':'k'})
# "distance" determines how far the numbers will be from the structure
# "shape" tunes the direction along which the numbers will move away from the structure
# "number spacing" can either be an integer or a list.
# An integer "i" will cause every "i-th" number to be displayed.
# A list will determine exactly what numbers to be shown.

# We can add a number at a distance=10, with shape=2 and spacing=4

nc.numbering_outside(10, 2, 4, {'fontsize':10, 'color':'k'})
```


    
![png](README_files/README_6_0.png)
    



```python
# Let's draw a 3-strands complex and color the strands differently

seq1 = 'UGACGUAAAACUGAC'
seq2 = 'UGUUACCGUAGUACG'
seq3 = 'ACCGUAC'
seq = "".join([seq1, seq2, seq3])
mfe = '..((((..(((....+.))))).)).(((((+..)))))'

nc = NucDraw(mfe)
nc.generate()
nc.plotter(8, bckwargs={'lw':2, 'color':'k'}, bpkwargs={'lw':3, 'c':'k'}, scwargs={'s':10, 'c':'k'})
nc.plot_circles(seq, circle_size = 2, circle_color='white')
nc.multistrand_coloring(clr=['red', 'blue', 'green'], bckwargs={'lw' : 3})

# Let's annotate with numbers in custom positions [3, 5, 8, 24].
# nc.numbering_outside(distance, shape, number spacing, {'fontsize':10, 'color':'k'})
# Notice that these positions may exceed the length of your strands, so message is thrown at the user.

nc.numbering_outside(10, 2, [3,5,8,24], {'fontsize':10, 'color':'k'})
```

    24 exceeds the length of the structure # 1 .
    8 exceeds the length of the structure # 2 .
    24 exceeds the length of the structure # 2 .



    
![png](README_files/README_7_1.png)
    



```python

```
