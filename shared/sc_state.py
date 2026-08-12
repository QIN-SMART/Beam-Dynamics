"""
SC runtime state contract (v0.14.1 task 3) — single source for both routes.

States:
  sc_requested   : the caller explicitly asked for SC physics (config
                   `space_charge.enabled` or an explicit function argument).
                   `step >= 4` only defines ROUTE CAPABILITY, not the
                   requested state.
  sc_available   : ocelot SpaceCharge import succeeded
  sc_configured  : SpaceCharge(...) constructed successfully
  sc_attached    : navi.add_physics_proc succeeded (coverage anchors set)
  sc_apply_count : number of apply() calls actually executed in tracking
  sc_effective   : = (sc_apply_count > 0)   ← the ONLY "SC actually ran" truth

HARD FAIL: when sc_requested is True, any of
    import fail / construct fail / attach fail / final apply_count == 0 /
    coverage != expected [cathode → sample]
raises immediately (RuntimeError).  There is NO silent fallback to no-SC:
a requested-but-not-running SC must never produce a normal-looking "SC ON"
result.

AG route: no scheduler states (available/configured/attached are OCELOT
concepts); AG metadata uses requested/effective + physical charge numbers.

Pure python (dataclasses only) — safe for shared/.
"""

from dataclasses import dataclass


@dataclass
class SCState:
    """OCELOT SC runtime state machine (one instance per run)."""
    requested: bool = False
    available: bool = False
    configured: bool = False
    attached: bool = False
    apply_count: int = 0
    coverage_start_m: float = 0.0
    coverage_stop_m: float = 0.0

    @property
    def effective(self) -> bool:
        """The only 'SC actually ran' truth."""
        return self.apply_count > 0

    def fail(self, msg: str):
        raise RuntimeError(f"SC runtime state HARD FAIL: {msg}")

    def check_attached(self):
        """HARD FAIL for import / construction / attach (call after setup)."""
        if not self.requested:
            return
        if not self.available:
            self.fail("sc_requested=True but SpaceCharge import failed "
                      "(unavailable)")
        if not self.configured:
            self.fail("sc_requested=True but SpaceCharge construction failed")
        if not self.attached:
            self.fail("sc_requested=True but add_physics_proc failed to "
                      "attach (no coverage anchors)")

    def verify_final(self, expected_start_m: float, expected_stop_m: float):
        """HARD FAIL after tracking: apply_count>0 and coverage == expected.

        expected_start_m / expected_stop_m : cathode → sample positions
        derived from lattice.elements (no hardcoded sample position).
        """
        if not self.requested:
            return
        self.check_attached()
        if self.apply_count <= 0:
            self.fail(f"sc_requested=True but sc_apply_count={self.apply_count} "
                      f"(==0) — SC did not run; refusing to emit a normal "
                      f"SC-ON result")
        if (abs(self.coverage_start_m - expected_start_m) > 1e-9 or
                abs(self.coverage_stop_m - expected_stop_m) > 1e-9):
            self.fail(f"sc_requested=True but coverage "
                      f"[{self.coverage_start_m}, {self.coverage_stop_m}] != "
                      f"expected [{expected_start_m}, {expected_stop_m}]")

    def to_meta(self) -> dict:
        """Flat metadata dict (recorded in the result meta)."""
        return {
            "sc_requested": self.requested,
            "sc_available": self.available,
            "sc_configured": self.configured,
            "sc_attached": self.attached,
            "sc_apply_count": self.apply_count,
            "sc_effective": self.effective,
            "sc_coverage_start_m": self.coverage_start_m,
            "sc_coverage_stop_m": self.coverage_stop_m,
        }
