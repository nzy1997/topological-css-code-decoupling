"""Shared finite-torus geometry."""

from dataclasses import dataclass


@dataclass(frozen=True)
class TorusShape:
    """Finite torus dimensions."""

    Lx: int
    Ly: int

    def __post_init__(self):
        """Validate that both torus periods are positive."""
        if self.Lx <= 0 or self.Ly <= 0:
            raise ValueError("Torus dimensions must be positive.")

    @property
    def size(self):
        """Return the number of torus sites."""
        return self.Lx * self.Ly

    def normalize(self, px, py):
        """Reduce a coordinate pair modulo the torus periods.

        Args:
            px: Finite torus x-coordinate.
            py: Finite torus y-coordinate.

        Returns:
            tuple: Normalized ``(px, py)`` coordinates.
        """
        return px % self.Lx, py % self.Ly

    def index(self, px, py):
        """Return the row-major vector index of a torus site.

        Args:
            px: Finite torus x-coordinate.
            py: Finite torus y-coordinate.

        Returns:
            int: Row-major site index.
        """
        ix, iy = self.normalize(px, py)
        return ix * self.Ly + iy

    def site(self, index):
        """Return ``(px, py)`` for one row-major site index.

        Args:
            index: Integer index into the relevant finite vector or table.

        Returns:
            tuple: ``(px, py)`` for the row-major site index.
        """
        return divmod(index, self.Ly)

    def sites(self):
        """Yield all ``(px, py)`` torus sites in row-major order."""
        for px in range(self.Lx):
            for py in range(self.Ly):
                yield px, py
