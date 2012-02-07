Constructing the methods
------------------------

First, start with creating the coefficients of the commutators in the BCH
expansion. This will come in handy in a minute. You actually do not need
to create them, the resultant files are available in this repository.
Here is how to create them though

     $ python RKNA_BCH_expan_coeffs.py
     $ python RKNB_BCH_expan_coeffs.py

This gives two files for your viewing 'D5_RKNA.dat' and 'D5_RKNB.dat',
and two more files for MAPLE's consumption 'RKN_methodA_coeffs.mw',
and 'RKN_methodB_coeffs.mw'.

Next, let's solve the order conditions. The equations are slightly
different for two types of RKN methods, arising from different splitting
arrangements. See the manuscript for details.

We first create MAPLE inputs with Python programs by running, e.g,

$ python random_grid_methA.py > random_grid_methA_seed17.mpl
$ python random_grid_methB.py > random_grid_methB_seed18.mpl

Modify random seeds (! and output file names ! both on command line and
in MAPLE programs created !) and number of tries appropriately.

Now we can feed these into MAPLE

$ maple < random_grid_methA_seed17.mpl
$ maple < random_grid_methB_seed18.mpl

This will produce some solutions and write them
into two files 'ord_cond_search_methA_seed17.mw' and
'ord_cond_search_methB_seed18.mw'. Incidentally, if you are not
happy with the MAPLE echoing back the contents of the script, replace
semicolons with colons in the Python code that creates mpl files.

The output is hard to parse for a human, so we reprint them in a more
readable form. At the same time, we increase the accuracy of the
solutions by setting the number of digits to 24. Finally, we also
calculate the error coefficients and classify the solutions. There are
three classes: PURELY IMAGINARY (leading order error terms are purely
imaginary), REAL (the solution to order conditions and hence the
timestep coefficients are purely real), OTHER. We do this in two steps
first a Python program creates a MAPLE program

$ python fancy_output_RKNA.py > fancy_output_RKNA_seed17.mpl
$ python fancy_output_RKNB.py > fancy_output_RKNB_seed18.mpl

These Python scripts have slightly different functionality. The second
one is capable of reading in multiple files. They also try to detect and
remove duplicate solutions. This functionality is not working perfectly
though, some duplicates remain.
We then feed the output into MAPLE:

$ maple < fancy_output_RKNA_seed17.mpl
$ maple < fancy_output_RKNB_seed18.mpl

Again, pay attention to changing random seeds at appropriate places. You
can easily overwrite a file that took a long time to compute. This is
bad programming on my part.

This should give you two files 'fancy_output_RKNA_seed17.dat' and
'fancy_output_RKNB_seed18.dat'. They are definitely human readable but
could still be a chore to sort through, especially since some duplicates
remain.

The output from my runs are in the repository with the names
'RKNA_solutions.dat' and 'RKNB_solutions.dat'. There may be more
solutions, especially in the RKNB family.

An Idea for Optimization
------------------------

See the manuscript for the idea on improving the obtained methods. I
used a similar methodology here, writing Python programs to create MAPLE
programs. I use the program 'trying_to_minimize_RKNB_6stage.py'. This
script takes values from already found solutions and feeds them into
MAPLE's minimization routine. I tried many different ways. Starting from
random values led to no solution at all. Taking linear combinations of
two RKNA solutions close to each other gave some promising results, but
the best results were obtained by starting in the vicinity of the best
solution, then taking the best solution from that search and starting in
its vicinity. This process can be largely automated. If you want to
monitor the outcome of the MAPLE optimization you can use something like
this

$ watch -n 12 'grep "^\[" minimization_search2_sixstage_RKNB_*.mw |cut -f 1 -d "," |cut -f2 -d "[" |sort -rg |tail '

Running this optimization requires BCH expansion
coefficients for a six-stage RKNB method. This is given in
'RKN_method_sixstageB_coeffs.mw' and produced by a simple modification
of the program mentioned at the beginning of this README. The result of
the "optimization" is in the file
'minimization_search_sixstage_RKNB_mock.mw'


Testing the methods
-------------------

Again, see the manuscript on the test methodology. The source
code provided here (in the directory testing) has been constantly
tinkered with while I was testing different methods. The version
here is with the scaffolding removed.  However, it still provides
all the functionality. In particular, the GBS integrator can be of
independent interest.

2-body, 400-body and timing runs are seperated into different
directories. There are shell (bash) scripts to make runs and Python
programs to read in the data and make plots.

For the timing run, you need 'err_check_mock'. You can get this
executable by modifying the 'Makefile', just change BNAME variable and
(re)run make.


Terms for Usage
---------------

This software is provided to help you get the job done. Hopefully
to solve a problem, but possibly as a nice pastime. There are
no restrictions whatsoever on the usage of this software. You can
modify and adapt it to any use you see fit. However, I cannot be held
responsible for any damage arising from using this software, including
but not limited to wasting your time. Also, distributing the software,
even parts of it as part of another software, is not usage. Such
actions are governed by "Terms for copying and distribution", as
indicated below.

I appreciate, but do not require, attribution. If you are unsure
about what reference to use for attribution, feel free to ask. If you
want to give attribution, proper reference would be Gürkan (2012,
submitted to Mathematics of Computation). Note that, if you simply
use this software to obtain some results to write a paper, you are
not required to offer me co-authorship or such. Only if I provide
extensive help to you for using the software, and hence I contribute
significant amount of my time specifically for your project, I would
like to be offered co-authorship. Do not let this to discourage you
from asking questions, but please read the documentation (at least
this page and the paper) before doing so.


Terms for Copying and Distribution
----------------------------------

This software is copyright (c) 2012 by Atakan Gürkan, except where
noted in the source files.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU [*General Public License version 3*](http://www.gnu.org/licenses/gpl-3.0.txt) as
published by the Free Software Foundation.

Note that I leave out the phrase "... or any later version", that is,
you cannot merge this software with software under different versions
of GPL. If you need to do this, please contact me and if the version
you propose is reasonable for me, I will make an addition or exception
to these terms.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

Parts of this software may be provided with licenses different from
but compatible with GNU General Public License version 3. All such
parts are clearly marked in the source code files. For those parts,
you are free to use GNU General Public License version 3 or the
particular license indicated in the source code.
