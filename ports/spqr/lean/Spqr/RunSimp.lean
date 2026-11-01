import Lean

/-- Simp set that unfolds `StateM` runs and the cost potentials of the walk / relabel phases. -/
register_simp_attr run_simp
register_simp_attr sim_simp
