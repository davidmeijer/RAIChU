# -*- coding: utf-8 -*-

"""This module provides RAIChU endpoints for RetroMol."""

import warnings

warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated",
    category=UserWarning
)

from typing import List, Optional, Union

from raichu.antismash import get_nrps_pks_modules
from raichu.cluster.modular_cluster import ModularCluster
from raichu.module import LinearPKSModule, IterativePKSModule, TransATPKSModule, NRPSModule
from raichu.representations import ClusterRepresentation
from raichu.run_raichu import build_cluster


BiosyntheticModule = Union[LinearPKSModule, IterativePKSModule, TransATPKSModule, NRPSModule]


def gbk_to_cluster_files(
    gbk_path: str, 
    out_path_core: str,
    out_path_tailoring: Optional[str] = None
) -> None:
    """
    Read functional modules and tailoring elements from GenBank cluster file.

    Make sure GenBank file contains a single cluster.
    
    :param gbk_path: Path to GenBank cluster file.
    :param out_path_core: Path to output core cluster file.
    :param out_path_tailoring: Path to output tailoring cluster file (optional).
    """
    module_blocks, cluster_info = get_nrps_pks_modules(gbk_path, return_info=True)
    cluster_repr = module_blocks.make_raichu_cluster()
    cluster_repr.write_cluster_core(out_path_core, cluster_info)
    if out_path_tailoring:
        cluster_repr.write_cluster_tailoring(out_path_tailoring)


def read_cluster_files(
    cluster_file_core: str,
    cluster_file_tailoring: Optional[str] = None
) -> Optional[ModularCluster]:
    """
    Parse cluster files and return a ModularCluster instance.

    :param cluster_file_core: Path to core cluster file.
    :param cluster_file_tailoring: Path to tailoring cluster file (optional).
    :return: A ModularCluster instance or None if cluster contains no functional modules.
    """
    cluster_repr = ClusterRepresentation.from_cluster_files(cluster_file_core, cluster_file_tailoring)
    cluster = build_cluster(cluster_repr, strict=False)
    
    # Check for functional modules
    has_functional_modules = False
    for module in cluster.modules:
        if not module.is_broken:
            has_functional_modules = True
    
    if not has_functional_modules:
        return None

    return cluster


def get_biosynthetic_modules(cluster: "ModularCluster") -> Optional[List["BiosyntheticModule"]]:
    """
    Return the longest contiguous sequence of modules such that:
      - The first module is a starter (is_starter_module == True) and NOT a terminator.
      - The last module is a terminator (is_termination_module == True) and NOT a starter.
      - Every module strictly between start and end is neither a starter nor a terminator.
    If no such sequence exists, return None.
    """
    modules = getattr(cluster, "modules", []) or []

    best_start = best_end = None  # inclusive indices of the best window
    current_start = None          # index of an active window's start (starter-only)
    # While a window is open, any flagged module (starter or terminator) can only appear as the valid end.
    # Otherwise it invalidates the window and (if starter-only) may start a new one.

    for i, m in enumerate(modules):
        is_start_only = bool(getattr(m, "is_starter_module", False)) and not bool(getattr(m, "is_termination_module", False))
        is_term_only  = bool(getattr(m, "is_termination_module", False)) and not bool(getattr(m, "is_starter_module", False))
        is_flagged    = bool(getattr(m, "is_starter_module", False) or getattr(m, "is_termination_module", False))

        if current_start is None:
            # Try to open a new window.
            if is_start_only:
                current_start = i
            # else remain closed
        else:
            # We have an open window; only neutral modules may appear in the middle.
            if is_term_only:
                # Valid close: compute length and update best if longer.
                if best_start is None or (i - current_start) > (best_end - best_start):
                    best_start, best_end = current_start, i
                current_start = None  # close window
            elif is_flagged:
                # Invalid middle (starter, terminator that's not term-only, or both flags True).
                # If this very module can start a new window, restart from here.
                current_start = i if is_start_only else None
            else:
                # Neutral module, continue window.
                pass

    if best_start is None:
        return None

    return modules[best_start:best_end + 1]


def draw_cluster(cluster: ModularCluster) -> str:
    """
    Draw cluster as SVG.

    :param cluster: The cluster to draw.
    :param out_path: The output path for the SVG file, 
    """
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()

    return cluster.draw_cluster(as_string=True)