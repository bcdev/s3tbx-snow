Some notes and instructions
------------------------------

'doc' is the s3tbx-snow documentation folder. Documentation is build from *.rst files using Sphinx.
index.rst is the main documentation page.

To install Sphinx run:

     > pip install Sphinx

To build the s3tbx-snow html documentation run:

     > cd ../s3tbx-snow/src/main/resources
     > sphinx-build -E -a -b html doc doc/_build/html

Then find the html documentation in ./doc/_build/html

Latex/pdf generation:
     > sphinx-build -E -a -b latex doc doc/_build/latex
     > cd doc/_build/latex
     > pdflatex SEOM_S3_SNOW_D2.3_SUM.tex

Then find the pdf file: ../s3tbx-snow/src/main/resources/doc/_build/latex/SEOM_S3_SNOW_D2.3_SUM.pdf

Change style elements (e.g. layout of titlepage):
- set 'attributes' in conf.py whereever possible
- set style definitions or new Latex commands in mystyle.sty
- redo Latex/pdf generation

Upload to readthedocs:
- just commit all s3tbx-snow changes to Github
- login at https://readthedocs.org/accounts/login/ with github user/passwd (dolaf/...)
- go to s3snow project main page
- if new build is not triggered automatically (should be?!), click button 'build version'. Then everything is done
automatically
- when new version is ready, go to downloads and check new pdf

More info:
    * Sphinx Tutorial: http://sphinx-doc.org/tutorial.html
    * RST Primer: http://sphinx-doc.org/rest.html#rst-primer