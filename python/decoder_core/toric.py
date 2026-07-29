"""Toric-code stacks decoded with PyMatching."""

from dataclasses import dataclass

from sage.all import GF, Matrix, vector


def _require_pymatching():
    """Import PyMatching with a clear error if it is unavailable.

    Returns:
        module: Imported ``pymatching`` module.
    """
    try:
        import pymatching
    except ImportError as exc:
        raise ValueError("pymatching is not installed in the active Sage Python environment.") from exc
    return pymatching


@dataclass(frozen=True)
class ToricMatchingDecoder:
    """Independent toric matching decoders for multiple stacks."""

    shape: object
    num_stacks: int
    weights: object = None

    def __post_init__(self):
        """Build one PyMatching decoder for each toric stack."""
        if self.num_stacks <= 0:
            raise ValueError("num_stacks must be positive.")
        pymatching = _require_pymatching()
        h = self._single_stack_check_matrix()
        h_rows = [[int(h[i, j]) for j in range(h.ncols())] for i in range(h.nrows())]
        matchings = tuple(
            pymatching.Matching(h_rows, weights=self.weights) for _ in range(self.num_stacks)
        )
        object.__setattr__(self, "_single_check_matrix", h)
        object.__setattr__(self, "_single_check_rows", h_rows)
        object.__setattr__(self, "_matchings", matchings)

    @property
    def vertices_per_stack(self):
        """Return the number of toric vertices per stack."""
        return self.shape.size

    @property
    def edges_per_stack(self):
        """Return the number of toric edges per stack."""
        return 2 * self.shape.size

    @property
    def syndrome_size(self):
        """Return the total syndrome length over all stacks."""
        return self.num_stacks * self.vertices_per_stack

    @property
    def correction_size(self):
        """Return the total edge-correction length over all stacks."""
        return self.num_stacks * self.edges_per_stack

    def vertex_index(self, stack, px, py):
        """Return total vertex index.

        Args:
            stack: Toric component stack index.
            px: Finite torus x-coordinate.
            py: Finite torus y-coordinate.

        Returns:
            int: Vertex index in the stacked finite layout.
        """
        return stack * self.vertices_per_stack + self.shape.index(px, py)

    def edge_index(self, stack, orientation, px, py):
        """Return total edge index.

        Args:
            stack: Toric component stack index.
            orientation: Toric edge orientation label.
            px: Finite torus x-coordinate.
            py: Finite torus y-coordinate.

        Returns:
            int: Edge index in the stacked finite layout.
        """
        if orientation not in ("h", "v"):
            raise ValueError("Toric edge orientation must be 'h' or 'v'.")
        offset = 0 if orientation == "h" else self.shape.size
        return stack * self.edges_per_stack + offset + self.shape.index(px, py)

    def edge_key(self, edge_index):
        """Return ``(stack, orientation, x, y)`` for a total edge index.

        Args:
            edge_index: Linear edge index in the stacked toric code.

        Returns:
            tuple: ``(stack, orientation, px, py)`` edge key.
        """
        stack, local = divmod(edge_index, self.edges_per_stack)
        if local < self.shape.size:
            orientation = "h"
            site_index = local
        else:
            orientation = "v"
            site_index = local - self.shape.size
        px, py = self.shape.site(site_index)
        return stack, orientation, px, py

    def zero_syndrome(self):
        """Return the zero stacked toric syndrome vector."""
        return vector(GF(2), self.syndrome_size)

    def zero_correction(self):
        """Return the zero stacked toric edge-correction vector."""
        return vector(GF(2), self.correction_size)

    def _single_stack_check_matrix(self):
        """Build one toric parity-check matrix."""
        # Horizontal edges connect (px, py) to (px + 1, py); vertical edges
        # connect (px, py) to (px, py + 1).  This ordering matches
        # ``edge_index`` and therefore the PyMatching weight vector order.
        h = Matrix(GF(2), self.vertices_per_stack, self.edges_per_stack, 0)
        for px, py in self.shape.sites():
            col_h = self.edge_index(0, "h", px, py)
            h[self.shape.index(px, py), col_h] += 1
            h[self.shape.index(px + 1, py), col_h] += 1

            col_v = self.edge_index(0, "v", px, py)
            h[self.shape.index(px, py), col_v] += 1
            h[self.shape.index(px, py + 1), col_v] += 1
        return h

    def _validate_length(self, values, expected, name):
        """Require an input vector to have the expected total length.

        Args:
            values: Input sequence to normalize or validate.
            expected: Expected vector length used for validation.
            name: Human-readable field name used in validation errors.
        """
        if len(values) != expected:
            raise ValueError(f"{name} length must be {expected}.")

    def boundary(self, correction):
        """Return the toric boundary of an edge correction.

        Args:
            correction: Candidate GF(2) correction vector.

        Returns:
            vector: The toric boundary of an edge correction.
        """
        self._validate_length(correction, self.correction_size, "correction")
        out = self.zero_syndrome()
        for stack in range(self.num_stacks):
            start_e = stack * self.edges_per_stack
            end_e = start_e + self.edges_per_stack
            local = vector(GF(2), correction[start_e:end_e])
            local_boundary = self._single_check_matrix * local
            start_v = stack * self.vertices_per_stack
            for i, value in enumerate(local_boundary):
                if value:
                    out[start_v + i] = 1
        return out

    def _validate_even_parity(self, syndrome):
        """Require each stack to have even syndrome parity.

        Args:
            syndrome: Stacked toric GF(2) syndrome vector.
        """
        for stack in range(self.num_stacks):
            start = stack * self.vertices_per_stack
            parity = (
                sum(
                    int(syndrome[start + i])
                    for i in range(self.vertices_per_stack)
                )
                % 2
            )
            if parity:
                raise ValueError("Each toric stack syndrome must have even parity.")

    def decode(self, syndrome):
        """Decode a multi-stack toric syndrome.

        Args:
            syndrome: Stacked toric GF(2) syndrome vector.

        Returns:
            vector: Stacked toric GF(2) edge correction.
        """
        self._validate_length(syndrome, self.syndrome_size, "syndrome")
        self._validate_even_parity(syndrome)
        out = self.zero_correction()
        for stack, matching in enumerate(self._matchings):
            start_v = stack * self.vertices_per_stack
            local_syndrome = [
                int(syndrome[start_v + i])
                for i in range(self.vertices_per_stack)
            ]
            local_correction = matching.decode(local_syndrome)
            start_e = stack * self.edges_per_stack
            for i, value in enumerate(local_correction):
                if value:
                    out[start_e + i] = 1
        return out

    def decode_with_weights(self, syndrome, weights):
        """Decode a multi-stack toric syndrome with per-edge weights.

        Args:
            syndrome: Stacked toric GF(2) syndrome vector.
            weights: Per-edge matching weights in stacked toric order.

        Returns:
            vector: Stacked toric GF(2) edge correction.
        """
        self._validate_length(syndrome, self.syndrome_size, "syndrome")
        self._validate_length(weights, self.correction_size, "weights")
        self._validate_even_parity(syndrome)
        pymatching = _require_pymatching()
        out = self.zero_correction()
        for stack in range(self.num_stacks):
            start_v = stack * self.vertices_per_stack
            local_syndrome = [
                int(syndrome[start_v + i])
                for i in range(self.vertices_per_stack)
            ]
            start_e = stack * self.edges_per_stack
            local_weights = [
                float(weights[start_e + i])
                for i in range(self.edges_per_stack)
            ]
            matching = pymatching.Matching(self._single_check_rows, weights=local_weights)
            local_correction = matching.decode(local_syndrome)
            for i, value in enumerate(local_correction):
                if value:
                    out[start_e + i] = 1
        return out
