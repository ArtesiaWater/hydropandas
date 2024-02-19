"""run notebooks in the examples directory."""

import os

import nbformat
import pytest
from nbconvert.preprocessors import ExecutePreprocessor

tst_dir = os.path.dirname(os.path.realpath(__file__))
nbdir = os.path.join(tst_dir, "..", "examples")


def _run_notebook(nbdir, fname):
    fname_nb = os.path.join(nbdir, fname)
    with open(fname_nb) as f:
        nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    out = ep.preprocess(nb, {"metadata": {"path": nbdir}})

    return out


@pytest.mark.notebooks
def test_run_notebook_01_groundwater_observations():
    _run_notebook(nbdir, "01_groundwater_observations.ipynb")


@pytest.mark.notebooks
def test_run_notebook_02_knmi_observations():
    _run_notebook(nbdir, "02_knmi_observations.ipynb")


@pytest.mark.notebooks
def test_run_notebook_03_hydropandas_and_pastas():
    _run_notebook(nbdir, "03_hydropandas_and_pastas.ipynb")


# not tested because notebook has intentional errors
# @pytest.mark.notebooks
# def test_run_notebook_04_merging_observations():
#     _run_notebook(nbdir, "04_merging_observations.ipynb")


@pytest.mark.notebooks
def test_run_notebook_05_bronhouderportaal_bro():
    _run_notebook(nbdir, "05_bronhouderportaal_bro.ipynb")
