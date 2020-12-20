Tests with rst syntax
=====================

Bulletted and numbered lists
++++++++++++++++++++++++++++
* This is a bulleted list.
* It has two items, the second
  item uses two lines.

1. This is a numbered list.
2. It has two items too.


1. This is a numbered list.
#. It has two items too.

Nested lists
++++++++++++

* this is
* a list

  * with a nested list
  * and some subitems

* and here the parent list continues

Definition lists
++++++++++++++++

term (up to a line of text)
   Definition of the term, which must be indented

   and can even consist of multiple paragraphs

next term
   Description.


Bla
+++

| These lines are
| broken exactly like in
| the source file.

Literal blocks
++++++++++++++

This is a normal text paragraph. The next paragraph is a code sample::

   It is not processed in any way, except
   that the indentation is removed.

   It can span multiple lines.

This is a normal text paragraph again.


Doctest blocks
++++++++++++++
>>> 1 + 1
2

Grid Tables
+++++++++++

+------------------------+------------+----------+----------+
| Header row, column 1   | Header 2   | Header 3 | Header 4 |
| (header rows optional) |            |          |          |
+========================+============+==========+==========+
| body row 1, column 1   | column 2   | column 3 | column 4 |
+------------------------+------------+----------+----------+
| body row 2             | ...        | ...      |          |
+------------------------+------------+----------+----------+


Simple Tables
+++++++++++++

=====  =====  =======
A      B      A and B
=====  =====  =======
False  False  False
True   False  False
False  True   False
True   True   True
=====  =====  =======


