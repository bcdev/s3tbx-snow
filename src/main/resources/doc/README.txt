'doc' is the CAWA documentation folder. Documentation is build from *.rst files using Sphinx.
index.rst is the main documentation page.

To install Sphinx run:

     > pip install Sphinx

To build the CAWA documentation run:

     > cd ../snap-cawa/src/main/resources/doc
     > make html

or to force regeneration of docs, run:

     > cd ../snap-cawa/src/main/resources
     > sphinx-build -E -a -b html doc doc/_build/html

Then find the html documentation in ../snap-cawa/src/main/resources/doc/_build/html

More info:
    * Sphinx Tutorial: http://sphinx-doc.org/tutorial.html
    * RST Primer: http://sphinx-doc.org/rest.html#rst-primer