[default]
all_states = state root mult sym dE_global_eV
by_mult = state_rel state root mult sym dE_gs_eV dE_gs_nm osc confdiffsw
docx = state_rel sym dE_gs_nm dE_gs_eV osc confdiff_strs weights

[nosym]
# When the molecule lacks any symmetry there is only one irrep (symmetry = 1)
# and we don't need any columns holding this value.
all_states = state mult dE_global_eV
by_mult = state_rel state root dE_gs_eV dE_gs_nm osc confdiffsw
docx = state_rel dE_gs_nm dE_gs_eV osc confdiff_strs weights

[nosym_occ]
# Use occupation strings like ..(π₁)²(π₂)²(π₃)²(π₄)²(d₁)ᵘ.. instead of
# differences between configurations for .docx-export.
all_states = state mult dE_global_eV
by_mult = state_rel state dE_gs_eV dE_gs_nm osc confdiffoccw_strs
docx = state_rel dE_gs_nm dE_gs_eV osc confdiffocc_strs weights
