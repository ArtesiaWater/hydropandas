==========
Contribute
==========

HydroPandas is an open source software project that depends on contributions
from the community. Everyone is welcome, each small contribution, no matter if
it is a fix of a typo in the documentation, bug report, an idea, or a question,
is valuable.

Questions, bug reports and feature requests
-------------------------------------------

If you have question about the use of HydroPandas, feel free to ask them in the
`GitHub Discussions <https://github.com/ArtesiaWater/hydropandas/discussions>`_.
Bugs and feature requests can be submitted via
`GitHub Issues <https://github.com/ArtesiaWater/hydropandas/issues>`_.

Version control, Git, and GitHub
--------------------------------

The code is hosted on GitHub. To contribute you will need to sign up for a free
GitHub account. We use Git for version control to allow many people to work
together on the project. If you have no experience with Git we recommend to
install `Github Desktop <https://desktop.github.com/>`_.

Contributing guidelines
-----------------------

Proposals for changes to the Hydropandas code base can be submitted via a pull
request. You can find a pull request (or PR) tutorial in the 
`GitHub's Help Docs. <https://help.github.com/articles/using-pull-requests/>`_.

There are roughly 6 steps for contributing to HydroPandas:

1. Fork the HydroPandas git repository
2. Create a development environment
3. Install HydroPandas dependencies
4. Make changes to code and add tests
5. Update the documentation
6. Submit a pull request

For pull request we use the following guidelines (similar to the 
`geopandas guidelines <https://geopandas.org/en/stable/community/contributing.html>`_):

- All existing tests should pass. Please make sure that the test suite passes,
both locally and on GitHub Actions. Status on Github Actions will be visible on
a pull request. To trigger a check, make a PR to your own fork.
- New functionality should include tests. Please write reasonable tests for your
code and make sure that they pass on your pull request.
- Classes, methods, functions, etc. should have docstrings. The first line of a
docstring should be a standalone summary. Parameters and return values should be
documented explicitly.
- Follow PEP 8 when possible. We use 
`Black <https://black.readthedocs.io/en/stable/>`_ and 
`Flake8 <http://flake8.pycqa.org/en/latest/>`_ to ensure a consistent code
format throughout the project.
- We use `isort <https://pycqa.github.io/isort/>` to automatically sort imports.
- We encourage backward compatability between HydroPandas versions but do not
ensure it (yet) because of the rapid changes to the code base.