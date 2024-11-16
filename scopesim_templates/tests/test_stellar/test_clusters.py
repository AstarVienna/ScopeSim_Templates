# -*- coding: utf-8 -*-
"""Tests for stellar.clusters"""

import pytest
from pytest import approx
from astropy import units as u
from astropy.table import Table

from spextra import Spextrum

from scopesim_templates.rc import Source
from scopesim_templates.stellar import cluster


@pytest.fixture(scope="module")
def basic_cluster():
    return cluster()


class TestCluster:
    @pytest.mark.usefixtures("basic_cluster")
    def test_it_works(self, basic_cluster):
        assert isinstance(basic_cluster, Source)

    @pytest.mark.usefixtures("basic_cluster")
    def test_spectra_are_correct_type(self, basic_cluster):
        assert all(isinstance(spec, Spextrum)
                   for spec in basic_cluster.spectra.values())

    @pytest.mark.usefixtures("basic_cluster")
    def test_spectra_waverange(self, basic_cluster):
        # TODO: don't know if it makes sense to test like this
        wvrng_min = [s.waverange.min().to(u.um).value
                     for s in basic_cluster.spectra.values()]
        wvrng_max = [s.waverange.max().to(u.um).value
                     for s in basic_cluster.spectra.values()]
        assert all(wave == approx(0.115) for wave in wvrng_min)
        assert all(wave == approx(2.5) for wave in wvrng_max)

    @pytest.mark.usefixtures("basic_cluster")
    def test_field_is_correct_type(self, basic_cluster):
        assert isinstance(basic_cluster.fields[0].field, Table)
        assert len(basic_cluster.fields[0]) > 0

    @pytest.mark.usefixtures("basic_cluster")
    def test_as_many_spectra_as_unique_sptypes(self, basic_cluster):
        uniques = set(basic_cluster.fields[0]["spec_types"])
        assert len(uniques) == len(basic_cluster.spectra)
