"""Results class."""


class Results(object):
    """Results objects saves all results from the graph."""

    def __init__(
        self,
        aromatic_carbons=None,
        charges=None,
        electrophilicities=None,
        local_electrophilicities=None,
        nucleophilicities=None,
        local_nucleophilicities=None,
        buried_volumina=None,
        feature_vectors=None,
        model_path=None,
        barriers=None,
        weights=None,
    ):
        self.aromatic_carbons = aromatic_carbons
        self.charges = charges
        self.electrophilicities = electrophilicities
        self.local_electrophilicities = local_electrophilicities
        self.nucleophilicities = nucleophilicities
        self.local_nucleophilicities = local_nucleophilicities
        self.buried_volumina = buried_volumina
        self.feature_vectors = feature_vectors
        self.model_path = model_path
        self.barriers = barriers
        self.weights = weights
