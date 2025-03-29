"""Test cluster."""

from astropy.table import Table
from synphot import SourceSpectrum

from scopesim_templates.micado import cluster
from scopesim_templates.rc import Source


class TestCluster:
    """Test cluster."""

    def test_cluster_returns_source_object(self):
        """Test cluster."""
        sky = cluster()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0].field, Table)
