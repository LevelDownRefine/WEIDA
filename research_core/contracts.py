"""Common algorithm output, without equating different scientific metrics."""
from dataclasses import dataclass
import numpy as np


@dataclass
class AlgorithmResult:
    method: str
    features: dict
    names: dict
    units: dict
    definitions: dict
    artifacts: dict
    parameters: dict

    def validate(self, n_scans):
        families = set(self.features)
        if not families or any(set(field) != families for field in
                               (self.names, self.units, self.definitions)):
            raise ValueError('Each feature family needs names, units and definition')
        for family, values in self.features.items():
            values = np.asarray(values, dtype=float)
            if (values.shape != (n_scans, len(self.names[family]))
                    or not len(self.names[family]) or np.isinf(values).any()
                    or len(set(self.names[family])) != len(self.names[family])):
                raise ValueError(f'Invalid scan x feature table: {family}')
            self.features[family] = values
        return self
